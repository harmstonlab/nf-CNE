import argparse
import glob
import gzip
import os
import subprocess
import sys
import bisect
from collections import defaultdict
from pathlib import Path
import re
import pdb

import bx.align.axt
import pandas as pd
import pyranges as pr
from tqdm import tqdm

parser = argparse.ArgumentParser()
parser.add_argument(
    "-i1",
    "--input_1",
    help="Pre-generated pairwise alignment file. Format should be .net.axt (.net.axt.gz also accepted)",
    required=False,
)
parser.add_argument(
    "-i2",
    "--input_2",
    help="Pre-generated pairwise alignment file. Format should be .net.axt (.net.axt.gz also accepted), should be inverse of input_1",
    required=False,
)
parser.add_argument(
    "-r1",
    "--reference_1",
    help="Path to either an indexed FASTA file or a directory containing an indexed FASTA file for the reference genome of input file 1",
    required=True,
)
parser.add_argument(
    "-r2",
    "--reference_2",
    help="Path to either an indexed FASTA file or a directory containing an indexed FASTA file for the reference genome of input file 2",
    required=True,
)
parser.add_argument(
    "-f1",
    "--filter_1",
    help="A bed file to use with respect to co-ordinates in input_1",
    required=False,
)
parser.add_argument(
    "-f2",
    "--filter_2",
    help="A bed file to use with respect to co-ordinates in input_2",
    required=False,
)
parser.add_argument(
    "-c",
    "--columns",
    type=int,
    help="The number of columns (bases) to include in a similarity window. Default=50",
    default=50,
)
parser.add_argument(
    "-id",
    "--identity",
    type=int,
    default=50,
    help="The identity threshold. At least this many bases within the column view must match to be counted. Default=50",
)

parser.add_argument(
    "-ci",
    "--chr_include",
    help="A regex that chromosomes must match to be included. Defaults to all cononical human, mouse, drosophila and zebrafish:\n^(chr[1-9][0-9]?|chr[XYM]|chr[2-3][LR]|chr4|super_[0-9]+|chrZ|chrW|chr[AB]|NC_[0-9]+)$",
    default="^(chr[1-9][0-9]?|chr[XYM]|chr[2-3][LR]|chr4|super_[0-9]+|chrZ|chrW|chr[AB]|NC_[0-9]+)$",
)

parser.add_argument(
    "-cx",
    "--chr_exclude",
    help="A regex that chromosomes must not match to be included. Defaults to exclude chrUn*, *_random, *_alt, *_fix, *_hap and NCBI scaffolds:\n(_random$|_alt$|_fix$|^chrUn|^GL|^KI|^JH|^KB|^NT|^NW|^NZ|_hap|_fix|_decoy)",
    default="(_random$|_alt$|_fix$|^chrUn|^GL|^KI|^JH|^KB|^NT|^NW|^NZ|_hap|_fix|_decoy)",
)


if len(sys.argv) <= 1:
    sys.exit(parser.print_help())

args = parser.parse_args()

input_file1 = args.input_1
input_file2 = args.input_2
reference1 = args.reference_1
reference2 = args.reference_2
filter1 = args.filter_1
filter2 = args.filter_2
columns = args.columns
identity = args.identity
INCLUDE_CHROM_RE = re.compile(args.chr_include)
EXCLUDE_CHROM_RE = re.compile(args.chr_exclude)

def chrom_exclude(chrom):
    return not (bool(INCLUDE_CHROM_RE.match(chrom)) and not EXCLUDE_CHROM_RE.search(chrom))

def load_bed_as_intervals(bedfile, seq_lengths):
    intervals_by_chr = {}
    rev_intervals_by_chr = {}
    with open(bedfile, "rt") as f:
        for line in f:
            if line.strip() == "":
                continue
            chrom, start, end, *rest = line.strip().split("\t")
            if chrom_exclude(chrom):
                continue

            start, end = int(start), int(end) - 1 

            if chrom not in intervals_by_chr:
                intervals_by_chr[chrom] = []
            intervals_by_chr[chrom].append((start, end))

            if chrom in seq_lengths:
                rev_start = seq_lengths[chrom] - end - 1
                rev_end = seq_lengths[chrom] - start - 1

                if chrom not in rev_intervals_by_chr:
                    rev_intervals_by_chr[chrom] = []
                rev_intervals_by_chr[chrom].append((rev_start, rev_end))

    for chrom in intervals_by_chr:
        intervals_by_chr[chrom].sort()
        merged = []
        current_start, current_end = intervals_by_chr[chrom][0]
        for s, e in intervals_by_chr[chrom][1:]:
            if s <= current_end + 1:
                current_end = max(current_end, e)
            else:
                merged.append((current_start, current_end))
                current_start, current_end = s, e
        merged.append((current_start, current_end))
        intervals_by_chr[chrom] = merged

    for chrom in rev_intervals_by_chr:
        rev_intervals_by_chr[chrom].sort()
        rev_merged = []
        rev_current_start, rev_current_end = rev_intervals_by_chr[chrom][0]
        for s, e in rev_intervals_by_chr[chrom][1:]:
            if s <= rev_current_end + 1:
                rev_current_end = max(rev_current_end, e)
            else:
                rev_merged.append((rev_current_start, rev_current_end))
                rev_current_start, rev_current_end = s, e
        rev_merged.append((rev_current_start, rev_current_end))
        rev_intervals_by_chr[chrom] = rev_merged

    return intervals_by_chr, rev_intervals_by_chr


def position_in_intervals(chrom, pos, intervals_by_chr):
    if chrom not in intervals_by_chr:
        return False
    intervals = intervals_by_chr[chrom]

    i = bisect.bisect_right(intervals, (pos, float("inf")))
    if i == 0:
        return False
    start, end = intervals[i - 1]

    return start <= pos <= end


if identity <= 0:
    sys.exit("Identity must be greater than 0")

if columns <= 0:
    sys.exit("Columns must be greater than 0")
if identity > columns:
    sys.exit("Identity must be less than or equal to columns")


if input_file1 == input_file2 and None not in (input_file1, input_file2):
    sys.exit("Both Input Files are the same. Different files must be provided")

if (input_file1 and not input_file2) or (input_file2 and not input_file1):
    sys.exit("Only one input file was provided. Two files must be provided")

if input_file1:
    if not os.path.exists(input_file1):
        sys.exit(f"Unable to find input file: {input_file1}")
    if not os.path.exists(input_file2):
        sys.exit(f"Unable to find input file: {input_file2}")
    elif not (input_file1.endswith(".net.axt") or input_file1.endswith(".net.axt.gz")):
        sys.exit(
            "Invalid extension for input file 1. Valid extensions are '.net.axt' and '.net.axt.gz'"
        )
    if input_file1.endswith(".net.axt.gz"):
        if not input_file2.endswith(".net.axt.gz"):
            sys.exit(
                "Input file 2 must be a .net.axt.gz file if input file 1 is a .net.axt.gz file"
            )
        open_func = gzip.open
    else:
        if input_file2.endswith(".net.axt.gz"):
            sys.exit(
                "Input file 2 must be a .net.axt file if input file 1 is a .net.axt file"
            )
        open_func = open

reference1 = (glob.glob(f"{reference1}/*.fai") + [reference1])[0]
reference2 = (glob.glob(f"{reference2}/*.fai") + [reference2])[0]

if not os.path.exists(reference1):
    sys.exit("Unable to locate index file for reference 1")

if not os.path.exists(reference2):
    sys.exit("Unable to locateindex file for reference 1")

seq_lengths_1 = {}
seq_lengths_2 = {}

with open(reference1, "rt") as f:
    for line in f:
        chrom, length, *_ = line.strip().split("\t")
        if chrom_exclude(chrom):
            continue
        seq_lengths_1[chrom] = int(length)

with open(reference2, "rt") as f:
    for line in f:
        chrom, length, *_ = line.strip().split("\t")
        if chrom_exclude(chrom):
            continue
        seq_lengths_2[chrom] = int(length)


if filter1:
    if not os.path.exists(filter1):
        sys.exit(f"Unable to find filter file: {filter1}")
    if filter1.endswith(".bed"):
        Filter_1, Rev_Filter_1 = load_bed_as_intervals(filter1, seq_lengths_1)
        Filter_1_name = "_FILTERED_" + os.path.basename(filter1).removesuffix(".bed")
    else:
        sys.exit("Invalid extension for filter file 1. Valid extension is '.bed'")
else:
    Filter_1, Rev_Filter_1 = dict(), dict()
    Filter_1_name = ""

if filter2:
    if not os.path.exists(filter2):
        sys.exit(f"Unable to find filter file: {filter2}")
    if filter2.endswith(".bed"):
        Filter_2, Rev_Filter_2 = load_bed_as_intervals(filter2, seq_lengths_2)
        Filter_2_name = "_FILTERED_" + os.path.basename(filter2).removesuffix(".bed")
    else:
        sys.exit("Invalid extension for filter file 2. Valid extension is '.bed'")
else:
    Filter_2, Rev_Filter_2 = dict(), dict()
    Filter_2_name = ""
    
bases = "ATCG"
# default dict as otherwise have to add every possible combination of N/n/- etc...
score_matrix = defaultdict(lambda: defaultdict(int))

for base in bases:
    for base1 in bases:
        if base == base1:
            score_matrix[base][base1] = 1
            score_matrix[base.lower()][base1] = 1
            score_matrix[base][base1.lower()] = 1
            score_matrix[base.lower()][base1.lower()] = 1

FINAL_outfile_1 = input_file1.replace(
    ".net.axt", f"_{identity}I_{columns}col{Filter_1_name}.bed"
).replace(".gz", "")
FINAL_outfile_2 = input_file2.replace(
    ".net.axt", f"_{identity}I_{columns}col{Filter_2_name}.bed"
).replace(".gz", "")

CEs_1vs2 = {}
CEs_1vs2["Chromosome"] = []
CEs_1vs2["Start"] = []
CEs_1vs2["End"] = []
CEs_1vs2["Name"] = []
CEs_1vs2["Score"] = []

CEs_2vs1 = {}
CEs_2vs1["Chromosome"] = []
CEs_2vs1["Start"] = []
CEs_2vs1["End"] = []
CEs_2vs1["Name"] = []
CEs_2vs1["Score"] = []

mapping_1 = defaultdict(str)
mapping_2 = defaultdict(str)
num_CEs = 0
# qStr = Query String = Sequence of query
# tStr = Template String = Sequence of template/reference
# tPosList = Template Position List = Array of what the genomic position is for each corresponding index in tStr (if tStr[i]=="-", tPosList[i]=="-1")
# asm1_chr = Assembly 1 Chromosome
# CE_start = How far into the record the Conserved Element Starts
# CE_end = How far into the record the Conserved Element Ends
# seq_lengths = dictionary of chr:sequence length of assembly2 in axt record as if strand=="-" position is chrom size minus position
# asm2_chr = Assembly 2 Chromosome
# asm2_starnd = Strand of Assembly 2 as if "-" need to flip co-ordinates
# 1PosList = Query Position List = Array of what the genomic position is for each corresponding index in qStr (if qStr[i]=="-", qPosList[i]=="-1")
# ROUND is just 1/2 if first/second file, as if first file it's asm1 vs asm2 but second is asm2 vs asm1 so need to swap template/query coordinates
def locate_and_format(
    qStr,
    tStr,
    tPosList,
    asm1_chr,
    CE_start,
    CE_end,
    seq_lengths,
    asm2_chr,
    asm2_strand,
    qPosList,
    ROUND,
):
    global num_CEs
    while CE_start not in qPosList.keys() or CE_end not in tPosList.keys():
        CE_start+=1
    while CE_end not in qPosList.keys() or CE_end not in tPosList.keys():
        CE_end-=1
    # trim no matching bases from start/end
    # redundant when identity==100%, but needed to tidy ends in all other cases
    while CE_start <= CE_end and score_matrix[qStr[CE_start]][tStr[CE_start]] <= 0:
        CE_start += 1
    while CE_end >= CE_start and score_matrix[qStr[CE_end]][tStr[CE_end]] <= 0:
        CE_end -= 1
    # Skip if trimming invalidated the interval
    if CE_start > CE_end:
        return

    if asm2_strand == "-":
        qSize = seq_lengths[asm2_chr]
        qStart = qSize - qPosList[CE_end] + 1
        qEnd = qSize - qPosList[CE_start] + 1
    else:
        qStart = qPosList[CE_start]
        qEnd = qPosList[CE_end]

    matching = sum(
        [
            t.upper() == q.upper()
            for t, q in zip(qStr[CE_start : CE_end + 1], tStr[CE_start : CE_end + 1])
        ]
    )
    length = CE_end - CE_start + 1
    score = round(matching / length * 100, 2)

    if ROUND == 1:
        CE_str_1 = f"{asm1_chr}:{tPosList[CE_start]}-{tPosList[CE_end]}"
        CE_str_2 = f"{asm2_chr}:{qStart}-{qEnd}"
        if CE_str_1 not in mapping_1:
                if CE_str_2 not in mapping_2:
                    num_CEs += 1
                    CE_number = f"CE_{num_CEs}"
                    mapping_1[CE_str_1] = CE_number
                    mapping_2[CE_str_2] = CE_number
                else:
                    CE_number = mapping_2[CE_str_2]
        else:
            CE_number = mapping_1[CE_str_1]

        CEs_1vs2["Chromosome"].append(asm1_chr)
        CEs_1vs2["Start"].append(tPosList[CE_start] - 1)  
        CEs_1vs2["End"].append(tPosList[CE_end])
        CEs_1vs2["Name"].append(CE_number)
        CEs_1vs2["Score"].append(score)

        CEs_2vs1["Chromosome"].append(asm2_chr)
        CEs_2vs1["Start"].append(qStart - 1)  
        CEs_2vs1["End"].append(qEnd)
        CEs_2vs1["Name"].append(CE_number)
        CEs_2vs1["Score"].append(score)

    elif ROUND == 2:
        CE_str_1 = f"{asm2_chr}:{qStart}-{qEnd}"
        CE_str_2 = f"{asm1_chr}:{tPosList[CE_start]}-{tPosList[CE_end]}"
        if CE_str_2 not in mapping_2:
            if CE_str_1 not in mapping_1:
                num_CEs += 1
                CE_number = f"CE_{num_CEs}"
                mapping_1[CE_str_1] = CE_number
                mapping_2[CE_str_2] = CE_number
            else:
                CE_number = mapping_1[CE_str_1]
        else:
            CE_number = mapping_2[CE_str_2]

        CEs_1vs2["Chromosome"].append(asm2_chr)
        CEs_1vs2["Start"].append(qStart - 1)  
        CEs_1vs2["End"].append(qEnd)
        CEs_1vs2["Name"].append(CE_number)
        CEs_1vs2["Score"].append(score)

        CEs_2vs1["Chromosome"].append(asm1_chr)
        CEs_2vs1["Start"].append(tPosList[CE_start] - 1)  
        CEs_2vs1["End"].append(tPosList[CE_end])
        CEs_2vs1["Name"].append(CE_number)
        CEs_2vs1["Score"].append(score)


# record is just bx.align.axt.Reader record object
# seq_lengths is dictionary of chr:sequence length of assembly2 in axt record as if strand=="-" position is chrom size minus position
# ROUND is just 1/2 if first/second file, as if first file it's asm1 vs asm2 but second is asm2 vs asm1 so need to swap template/query coordinates
def scanAxt(record, seq_lengths, ROUND, tFilter, qFilter):
    asm1_chr = record.components[0].src.split(".")[-1]
    asm2_chr = record.components[1].src.split(".")[-1]
    asm2_strand = record.components[1].strand

    if (chrom_exclude(asm1_chr)) or (chrom_exclude(asm2_chr)):
        return None

    tStr, qStr = zip(*record.column_iter())
    tPos = record.components[0].start
    qPos = record.components[1].start

    profile = defaultdict(int)
    tPosList = {}
    qPosList = {}
    CE_start = None
    CE_end = None
    i = 0

    while i < len(tStr):
        # Skip over target filter intervals
        while position_in_intervals(asm1_chr, tPos, tFilter) or position_in_intervals(asm2_chr, qPos, qFilter):
            if tStr[i] != "-":
                tPos += 1
            if qStr[i] != "-":
                qPos += 1
            i += 1
            if i >= len(tStr):
                if CE_start is not None:
                    locate_and_format(
                        qStr, tStr, tPosList,
                        asm1_chr, CE_start, CE_end,
                        seq_lengths, asm2_chr, asm2_strand,
                        qPosList, ROUND
                    )
                    CE_start = CE_end = None
                return

        profile[i] = score_matrix[tStr[i]][qStr[i]]
        tPosList[i] = -1 if tStr[i] == "-" else (tPos := tPos + 1)
        qPosList[i] = -1 if qStr[i] == "-" else (qPos := qPos + 1)
        columns_seen = 1
        i += 1

        while i < len(tStr):
            while position_in_intervals(asm1_chr, tPos, tFilter) or position_in_intervals(asm2_chr, qPos, qFilter):
                if tStr[i] != "-":
                    tPos += 1
                if qStr[i] != "-":
                    qPos += 1
                i += 1
                if i >= len(tStr):
                    if CE_start is not None:
                        locate_and_format(
                            qStr, tStr, tPosList,
                            asm1_chr, CE_start, CE_end,
                            seq_lengths, asm2_chr, asm2_strand,
                            qPosList, ROUND
                        )
                        CE_start = CE_end = None
                    return

            tPosList[i] = -1 if tStr[i] == "-" else (tPos := tPos + 1)
            qPosList[i] = -1 if qStr[i] == "-" else (qPos := qPos + 1)
            profile[i] = profile[i - 1] + score_matrix[tStr[i]][qStr[i]]
            columns_seen += 1

            if columns_seen >= columns:
                score = profile[i] if columns_seen == columns else profile[i] - profile[i - columns]

                if score >= identity:
                    if CE_start is None:
                        CE_start = i - columns + 1
                    CE_end = i
                else:
                    if CE_start is not None and CE_end < i - columns + 1:
                        # Finalize and output CE
                        locate_and_format(
                            qStr, tStr, tPosList,
                            asm1_chr, CE_start, CE_end,
                            seq_lengths, asm2_chr, asm2_strand,
                            qPosList, ROUND
                        )
                        CE_start = CE_end = None
            i += 1

        if CE_start is not None:
            locate_and_format(
                qStr, tStr, tPosList,
                asm1_chr, CE_start, CE_end,
                seq_lengths, asm2_chr, asm2_strand,
                qPosList, ROUND
            )
            CE_start = CE_end = None
    if CE_start is not None:
        locate_and_format(
            qStr, tStr, tPosList,
            asm1_chr, CE_start, CE_end,
            seq_lengths, asm2_chr, asm2_strand,
            qPosList, ROUND
        )
        CE_start = CE_end = None





print(f"{f'Processing {os.path.basename(input_file1)}':-^60}")

with open_func(input_file1, "rt") as f:
    try:
        last_line = subprocess.check_output(["tail", "-4", input_file1])
        last_record = int(last_line.decode().split("\n")[0].split(" ")[0])
    except UnicodeDecodeError:
        lines = f.read().splitlines()
        last_record = int(lines[-4].split("\n")[0].split(" ")[0])
        f.seek(0)

    try:
        axt_data = bx.align.axt.Reader(f)
        for record in tqdm(axt_data, total=last_record, leave=False):
            if record.components[0].strand == "+":
                tFilter_to_use = Filter_1
            else:
                tFilter_to_use = Rev_Filter_1
            if record.components[1].strand == "+":
                qFilter_to_use = Filter_2
            else:
                qFilter_to_use = Rev_Filter_2
            scanAxt(
                record,
                seq_lengths_2,
                ROUND=1,
                tFilter=tFilter_to_use,
                qFilter=qFilter_to_use,
            )

    except ValueError as e:
        print(e)
        sys.exit("Unable to read provided axt file")

print(f"{f'Processing {os.path.basename(input_file2)}':-^60}")
with open_func(input_file2, "rt") as f:
    try:
        last_line = subprocess.check_output(["tail", "-4", input_file2])
        last_record = int(last_line.decode().split("\n")[0].split(" ")[0])
    except UnicodeDecodeError:
        lines = f.read().splitlines()
        last_record = int(lines[-4].split("\n")[0].split(" ")[0])
        f.seek(0)

    try:
        axt_data = bx.align.axt.Reader(f)
        for record in tqdm(axt_data, total=last_record, leave=False):
            if record.components[0].strand == "+":
                tFilter_to_use = Filter_2
            else:
                tFilter_to_use = Rev_Filter_2
            if record.components[1].strand == "+":
                qFilter_to_use = Filter_1
            else:
                qFilter_to_use = Rev_Filter_1
            scanAxt(
                record,
                seq_lengths_1,
                ROUND=2,
                tFilter=tFilter_to_use,
                qFilter=qFilter_to_use,
            )

    except ValueError:
        sys.exit("Unable to read provided axt file")

df_1vs2 = pd.DataFrame(CEs_1vs2)
df_2vs1 = pd.DataFrame(CEs_2vs1)
py_range_1vs2=pr.PyRanges(df_1vs2)
py_range_2vs1=pr.PyRanges(df_2vs1)

df_1vs2 = df_1vs2.rename(columns={
    "Chromosome": "Chromosome_species_1",
    "Start": "Start_species_1",
    "End": "End_species_1",
    "Score": "Score_species_1"
})

df_2vs1 = df_2vs1.rename(columns={
    "Chromosome": "Chromosome_species_2",
    "Start": "Start_species_2",
    "End": "End_species_2",
    "Score": "Score_species_2"
})

merged_df = pd.merge(df_1vs2, df_2vs1, on="Name", suffixes=('_1', '_2'))
df_species_1_coords = df_1vs2.rename(columns={
    "Chromosome_species_1": "Chromosome",
    "Start_species_1": "Start",
    "End_species_1": "End"
})
py_range_species_1 = pr.PyRanges(df_species_1_coords)

# Self-join to find contained elements in species_1
overlaps = py_range_species_1.join(py_range_species_1, suffix="_b")
if overlaps:
    # Convert to DataFrame and filter for strict containment in species_1
    overlap_df = overlaps.df
    contained_in_species_1 = overlap_df[
        (overlap_df.Start_b <= overlap_df.Start) &
        (overlap_df.End_b >= overlap_df.End) &
        ~((overlap_df.Start == overlap_df.Start_b) & (overlap_df.End == overlap_df.End_b))
    ]

    # Join with merged_df to get species_2 coordinates for both Name and Name_b
    contained = contained_in_species_1.merge(
        merged_df, on="Name"
    ).merge(
        merged_df[["Name", "Chromosome_species_2", "Start_species_2", "End_species_2"]],
        left_on="Name_b", right_on="Name", suffixes=('', '_b')
    )

    # Filter for cases where CNE is also nested in species_2 coordinates
    contained_in_both_species = contained[
        (contained.Chromosome_species_2 == contained.Chromosome_species_2_b) &
        (contained.Start_species_2 >= contained.Start_species_2_b) &
        (contained.End_species_2 <= contained.End_species_2_b)
    ]

    # Remove these contained names from df_1vs2
    contained_names = set(contained_in_both_species["Name"])
    filtered_df = df_1vs2[~df_1vs2["Name"].isin(contained_names)]

    filtered_df_for_pyranges = filtered_df.rename(columns={
        "Chromosome_species_1": "Chromosome",
        "Start_species_1": "Start",
        "End_species_1": "End"
    })
    filtered_py_range = pr.PyRanges(filtered_df_for_pyranges)

    # Get the list of kept names
    kept_names = set(filtered_df["Name"])

    # Filter merged_df using kept names
    merged_filtered = merged_df[merged_df["Name"].isin(kept_names)]

    # Extract species_1 BED columns
    bed_1vs2 = merged_filtered[[
        "Chromosome_species_1", "Start_species_1", "End_species_1", "Name", "Score_species_1"
    ]]
    bed_1vs2 = bed_1vs2.sort_values(["Chromosome_species_1", "Start_species_1"])

    # Extract species_2 BED columns
    bed_2vs1 = merged_filtered[[
        "Chromosome_species_2", "Start_species_2", "End_species_2", "Name", "Score_species_2"
    ]]
    bed_2vs1 = bed_2vs1.sort_values(["Chromosome_species_2", "Start_species_2"])

    # Write to BED files
    panda_CEs = pd.DataFrame(bed_1vs2.rename(columns={
        "Chromosome_species_1": "Chromosome",
        "Start_species_1": "Start",
        "End_species_1": "End",
        "Score_species_1":"Score"
    }))
    py_range_data = pr.PyRanges(panda_CEs)
    uniques = py_range_data.drop_duplicate_positions()
    print(f"{f'Saving output as {os.path.basename(FINAL_outfile_1)}':-^60}")
    if os.path.exists(FINAL_outfile_1):
        os.unlink(FINAL_outfile_1)
    if uniques:
        uniques.to_bed(FINAL_outfile_1)
    else:
        Path(FINAL_outfile_1).touch()

    panda_CEs = pd.DataFrame(bed_2vs1.rename(columns={
        "Chromosome_species_2": "Chromosome",
        "Start_species_2": "Start",
        "End_species_2": "End",
        "Score_species_2":"Score"
    }))
    py_range_data = pr.PyRanges(panda_CEs)
    uniques = py_range_data.drop_duplicate_positions()
    print(f"{f'Saving output as {os.path.basename(FINAL_outfile_2)}':-^60}")
    if os.path.exists(FINAL_outfile_2):
        os.unlink(FINAL_outfile_2)
    if uniques:
        uniques.to_bed(FINAL_outfile_2)
    else:
        Path(FINAL_outfile_2).touch()
else:
    Path(FINAL_outfile_1).touch()
    Path(FINAL_outfile_2).touch()
