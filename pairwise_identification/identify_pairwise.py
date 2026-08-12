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

import bx.align.axt
import pandas as pd
import pyranges as pr
from tqdm import tqdm

parser = argparse.ArgumentParser()
parser.add_argument(
    "-i1",
    "--input_1",
    help="Pre-generated pairwise alignment file. Format should be .net.axt (.net.axt.gz also accepted)",
    required=True,
)
parser.add_argument(
    "-i2",
    "--input_2",
    help="Pre-generated pairwise alignment file. Format should be .net.axt (.net.axt.gz also accepted), should be inverse of input_1",
    required=True,
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
    help="A regex that chromosomes must match to be included. Defaults to all canonical human, mouse, drosophila and zebrafish:\n^(chr[1-9][0-9]?|chr[XYM]|chr[2-3][LR]|chr4|chrLG[0-9]+|super_[0-9]+|chrZ|chrW|chr[AB]|NC_[0-9]+)$",
    default="^(chr[1-9][0-9]?|chr[XYM]|chr[2-3][LR]|chr4|chrLG[0-9]+|super_[0-9]+|chrZ|chrW|chr[AB]|NC_[0-9]+)$",
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


# input_1/input_2 are required by argparse, so both are always present here.
if input_file1 == input_file2:
    sys.exit("Both Input Files are the same. Different files must be provided")

if not os.path.exists(input_file1):
    sys.exit(f"Unable to find input file: {input_file1}")
if not os.path.exists(input_file2):
    sys.exit(f"Unable to find input file: {input_file2}")

if not (input_file1.endswith(".net.axt") or input_file1.endswith(".net.axt.gz")):
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
if not reference1.endswith(".fai"):
    sys.exit(f"Reference 1 is not a .fai index file: {reference1}")

if not os.path.exists(reference2):
    sys.exit("Unable to locate index file for reference 2")
if not reference2.endswith(".fai"):
    sys.exit(f"Reference 2 is not a .fai index file: {reference2}")

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
CEs_1vs2["Strand"] = []

CEs_2vs1 = {}
CEs_2vs1["Chromosome"] = []
CEs_2vs1["Start"] = []
CEs_2vs1["End"] = []
CEs_2vs1["Name"] = []
CEs_2vs1["Score"] = []
CEs_2vs1["Strand"] = []

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
    # Nudge the boundaries onto columns that were actually recorded (safety net now
    # that filtered columns break a segment; each loop tests the variable it moves).
    while CE_start <= CE_end and (CE_start not in qPosList or CE_start not in tPosList):
        CE_start += 1
    while CE_end >= CE_start and (CE_end not in qPosList or CE_end not in tPosList):
        CE_end -= 1
    if CE_start > CE_end:
        return
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
        q.upper() == t.upper()
        for q, t in zip(qStr[CE_start : CE_end + 1], tStr[CE_start : CE_end + 1])
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
        # Provisional per-pair gauge: anchor species_1 (A) on "+" so species_2 (B) can
        # carry the pair's relative orientation (the axt query strand). Only valid for
        # this row alone -- assign_consistent_strands rewrites both per locus later.
        CEs_1vs2["Strand"].append("+")

        CEs_2vs1["Chromosome"].append(asm2_chr)
        CEs_2vs1["Start"].append(qStart - 1)
        CEs_2vs1["End"].append(qEnd)
        CEs_2vs1["Name"].append(CE_number)
        CEs_2vs1["Score"].append(score)
        CEs_2vs1["Strand"].append(asm2_strand)

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
        # species_1 (A) is the query here but keeps the same provisional gauge; the
        # relative orientation is the same physical quantity, so species_2 (B) again
        # gets asm2_strand. Rewritten per locus by assign_consistent_strands.
        CEs_1vs2["Strand"].append("+")

        CEs_2vs1["Chromosome"].append(asm1_chr)
        CEs_2vs1["Start"].append(tPosList[CE_start] - 1)
        CEs_2vs1["End"].append(tPosList[CE_end])
        CEs_2vs1["Name"].append(CE_number)
        CEs_2vs1["Score"].append(score)
        CEs_2vs1["Strand"].append(asm2_strand)


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
    # Number of contiguous *kept* columns since the last filtered break. A filtered
    # column resets this to 0, so a similarity window never spans filtered bases and
    # the index arithmetic below (profile[i - columns], i - columns + 1) only ever
    # references columns from the current unbroken segment.
    columns_seen = 0

    def finalize():
        nonlocal CE_start, CE_end
        if CE_start is not None:
            locate_and_format(
                qStr, tStr, tPosList,
                asm1_chr, CE_start, CE_end,
                seq_lengths, asm2_chr, asm2_strand,
                qPosList, ROUND
            )
            CE_start = CE_end = None

    i = 0
    while i < len(tStr):
        # A filtered column is a hard break: finalize any open CE and reset the
        # segment so the next window starts fresh past the filtered region.
        if position_in_intervals(asm1_chr, tPos, tFilter) or position_in_intervals(asm2_chr, qPos, qFilter):
            finalize()
            columns_seen = 0
            if tStr[i] != "-":
                tPos += 1
            if qStr[i] != "-":
                qPos += 1
            i += 1
            continue

        # Record the genomic position of this column (1-based; -1 marks a gap).
        tPosList[i] = -1 if tStr[i] == "-" else (tPos := tPos + 1)
        qPosList[i] = -1 if qStr[i] == "-" else (qPos := qPos + 1)

        # Cumulative match score within the current segment. columns_seen == 0 means
        # this is the segment's first column, so seed the baseline non-cumulatively.
        if columns_seen == 0:
            profile[i] = score_matrix[tStr[i]][qStr[i]]
        else:
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
                    finalize()
        i += 1

    finalize()





def count_axt_records(path):
    """Best-effort total record count for the tqdm bar.

    The last axt record number lives on the summary line 4 lines from EOF. `tail`
    only works on the plaintext .net.axt; on a gzipped file it reads compressed
    bytes, so any failure (including a garbled decode) just yields None and the
    progress bar runs without a total.
    """
    try:
        last_line = subprocess.check_output(["tail", "-4", path])
        return int(last_line.decode().split("\n")[0].split(" ")[0])
    except UnicodeDecodeError:
        lines = f.read().splitlines()
        last_record = int(lines[-4].split("\n")[0].split(" ")[0])
        f.seek(0)
        return last_record
    except:
        return None


def assign_consistent_strands(merged_df):
    """Replace the per-pair strand gauge with one strand per LOCUS (edits in place).

    An axt record only carries a pair's RELATIVE orientation, so locate_and_format
    encodes it with a fixed gauge: species_1 written "+", species_2 carrying the axt
    query strand. Unambiguous for one row alone, but it breaks the moment a locus is
    reached by several pairs -- mm10 chr12:48289297-48291823 aligns to hg38 chr14
    (r=+), chr2 (r=-) and chr9 (r=+). With species_1 pinned to "+" the chr2 inversion
    is unrepresentable and is lost at write time.

    So treat each row as an edge of a bipartite graph (species_1 loci vs species_2
    loci) labelled with that relative orientation, and solve for one strand per locus
    satisfying strand_b = strand_a * r on every edge. Extracting every locus of a
    connected group on its listed strand then yields the same sequence.

    Each component's gauge is fixed by its highest-identity row (ties broken by append
    order, -i1 before -i2) -- the same priority the per-locus collapse further down
    uses -- with that row's species_1 locus on "+". Rows contradicting an
    already-fixed gauge are counted and ignored rather than overwriting the
    higher-identity assignment.

    Returns the number of contradicting rows.
    """
    if merged_df.empty:
        return 0

    # One integer id per distinct locus. groupby().ngroup() deduplicates in C, and
    # integer ids let the union-find state live in flat lists rather than dicts keyed
    # on millions of (chrom, start, end) tuples. dropna=False so a malformed row can
    # never land in the -1 "no group" bucket and silently share a node.
    a_id = merged_df.groupby(
        ["Chromosome_species_1", "Start_species_1", "End_species_1"],
        sort=False, dropna=False,
    ).ngroup().to_numpy()
    b_id = merged_df.groupby(
        ["Chromosome_species_2", "Start_species_2", "End_species_2"],
        sort=False, dropna=False,
    ).ngroup().to_numpy()

    # Both species share one id space; the offset keeps them distinct even when a
    # species_1 and a species_2 locus happen to share chromosome name and coordinates.
    n_a = int(a_id.max()) + 1
    n_nodes = n_a + int(b_id.max()) + 1

    # Union-find with parity. parity[x] is x's orientation relative to parent[x]:
    # 0 = same strand, 1 = reverse complement.
    parent = list(range(n_nodes))
    parity = bytearray(n_nodes)
    # Priority rank of the row that first created each root. Merges always hang the
    # junior root under the senior one, so a component is never re-rooted and its
    # anchor keeps "+" however many pairs join later. Union-by-size would re-root and
    # silently flip an already-anchored component.
    anchor = [-1] * n_nodes

    def find(x):
        """(root, parity of x relative to root), compressing the path walked."""
        root, p = x, 0
        while parent[root] != root:
            p ^= parity[root]
            root = parent[root]
        cur, cur_p = x, p
        while cur != root:
            nxt, nxt_p = parent[cur], cur_p ^ parity[cur]
            parent[cur], parity[cur] = root, cur_p
            cur, cur_p = nxt, nxt_p
        return root, p

    # Score_species_1 == Score_species_2 by construction (locate_and_format appends the
    # same identity to both), so either orders the rows identically. Keep it to a
    # SINGLE sort key: pandas implements descending sort as reverse -> stable sort ->
    # reverse, so equal scores keep their original append order (-i1 rows first), but
    # multiple keys route through lexsort_indexer instead and lose that property.
    edges = pd.DataFrame({
        "a": a_id,
        "b": b_id + n_a,
        "flip": merged_df["Strand_species_2"].to_numpy() == "-",
        "score": merged_df["Score_species_1"].to_numpy(),
    }).sort_values("score", ascending=False, kind="stable")

    conflicts = 0
    for rank, edge in enumerate(edges.itertuples(index=False)):
        a, b, flip = int(edge.a), int(edge.b), int(edge.flip)
        if anchor[a] < 0:
            anchor[a] = rank
        if anchor[b] < 0:
            anchor[b] = rank
        ra, pa = find(a)
        rb, pb = find(b)
        if ra == rb:
            # A cycle: only a problem if it demands the opposite of what is already set.
            if (pa ^ pb) != flip:
                conflicts += 1
        elif anchor[ra] <= anchor[rb]:
            parent[rb], parity[rb] = ra, pa ^ pb ^ flip
        else:
            parent[ra], parity[ra] = rb, pa ^ pb ^ flip

    # Resolve each locus once, then fan the per-locus strand back out to its rows.
    # Assigning a raw ndarray to a column is positional, so this is index-safe.
    lut = pd.Series(["-" if find(n)[1] else "+" for n in range(n_nodes)])
    merged_df["Strand_species_1"] = lut.take(a_id).to_numpy()
    merged_df["Strand_species_2"] = lut.take(b_id + n_a).to_numpy()
    return conflicts


print(f"{f'Processing {os.path.basename(input_file1)}':-^60}")

with open_func(input_file1, "rt") as f:
    last_record = count_axt_records(input_file1)

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
    last_record = count_axt_records(input_file2)

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

df_1vs2 = df_1vs2.rename(columns={
    "Chromosome": "Chromosome_species_1",
    "Start": "Start_species_1",
    "End": "End_species_1",
    "Score": "Score_species_1",
    "Strand": "Strand_species_1"
})

df_2vs1 = df_2vs1.rename(columns={
    "Chromosome": "Chromosome_species_2",
    "Start": "Start_species_2",
    "End": "End_species_2",
    "Score": "Score_species_2",
    "Strand": "Strand_species_2"
})

# CEs_1vs2 and CEs_2vs1 are filled in lockstep by locate_and_format: row i of one is
# the species-2 partner of row i of the other. Pair them by row position instead of
# pd.merge(on="Name"): a CE_n can be reused for several genuinely distinct element
# pairs (real multi-copy / paralogous hits), so a name-merge is many-to-many and would
# Cartesian-explode and mispair unrelated species_1 / species_2 coordinates.
merged_df = pd.concat(
    [df_1vs2.reset_index(drop=True), df_2vs1.drop(columns="Name").reset_index(drop=True)],
    axis=1,
)
# The per-pair "+"/relative-orientation gauge written by locate_and_format is only
# valid for a row in isolation; rewrite it as one strand per locus so every locus
# reachable through a chain of alignments agrees. Must run here -- before the
# containment self-join and before any filtering -- because a row dropped later can be
# the only bridge in a component, and gauging the halves separately would let them
# disagree.
strand_conflicts = assign_consistent_strands(merged_df)
if strand_conflicts:
    print(
        f"WARNING: {strand_conflicts} alignment pair(s) contradicted an already-fixed "
        "strand for a shared locus. The input alignments are mutually inconsistent "
        "there; kept the higher-identity orientation."
    )
# Build the self-join input from merged_df (not df_1vs2) so each species_1 interval
# carries its own row-paired species_2 partner coordinates as metadata. This lets the
# self-join below attach both partners (unsuffixed = contained, _b = container) without
# a merge on "Name" -- a name is reused across paralogous pairs, so merging on it is
# many-to-many and Cartesian-explodes (see the note above merged_df).
df_species_1_coords = merged_df.rename(columns={
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

    # The species_2 partner coordinates for both the contained element (unsuffixed) and
    # the container (_b) travelled through the self-join above, so filter for cases where
    # the CNE is also nested in species_2 coordinates directly -- no merge on "Name" needed.
    contained_in_both_species = contained_in_species_1[
        (contained_in_species_1.Chromosome_species_2 == contained_in_species_1.Chromosome_species_2_b) &
        (contained_in_species_1.Start_species_2 >= contained_in_species_1.Start_species_2_b) &
        (contained_in_species_1.End_species_2 <= contained_in_species_1.End_species_2_b)
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
        "Chromosome_species_1", "Start_species_1", "End_species_1", "Name", "Score_species_1", "Strand_species_1"
    ]].copy()
    # One locus can be reached by several alignment pairs. Strand is no longer at stake here
    # (assign_consistent_strands already gave every row sharing a locus the same strand), so this
    # collapse decides which Score and Name reach the BED. Keep the highest-identity pair; on a
    # tie keep the first one encountered (append order == file order, -i1 before -i2). pandas
    # implements descending sort as reverse -> stable sort -> reverse, so equal scores keep their
    # original append order.
    bed_1vs2 = bed_1vs2.sort_values("Score_species_1", ascending=False, kind="stable")
    bed_1vs2 = bed_1vs2.drop_duplicates(
        subset=["Chromosome_species_1", "Start_species_1", "End_species_1"], keep="first"
    )
    bed_1vs2 = bed_1vs2.sort_values(["Chromosome_species_1", "Start_species_1"])

    # Extract species_2 BED columns
    bed_2vs1 = merged_filtered[[
        "Chromosome_species_2", "Start_species_2", "End_species_2", "Name", "Score_species_2", "Strand_species_2"
    ]].copy()
    # Same deterministic per-locus collapse as species_1 above (highest identity, then first).
    bed_2vs1 = bed_2vs1.sort_values("Score_species_2", ascending=False, kind="stable")
    bed_2vs1 = bed_2vs1.drop_duplicates(
        subset=["Chromosome_species_2", "Start_species_2", "End_species_2"], keep="first"
    )
    bed_2vs1 = bed_2vs1.sort_values(["Chromosome_species_2", "Start_species_2"])

    # Write to BED files
    panda_CEs = pd.DataFrame(bed_1vs2.rename(columns={
        "Chromosome_species_1": "Chromosome",
        "Start_species_1": "Start",
        "End_species_1": "End",
        "Score_species_1":"Score",
        "Strand_species_1": "Strand"
    }))
    py_range_data = pr.PyRanges(panda_CEs)
    uniques = py_range_data.drop_duplicate_positions(strand=False)
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
        "Score_species_2":"Score",
        "Strand_species_2": "Strand"
    }))
    py_range_data = pr.PyRanges(panda_CEs)
    uniques = py_range_data.drop_duplicate_positions(strand=False)
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
