import argparse
import os
import sys
from collections import defaultdict

parser = argparse.ArgumentParser()
parser.add_argument("-b", "--bed", help="Bed file", required=True)
parser.add_argument(
    "-c",
    "--cutoff",
    help="Maximum number off hits allowed",
    default=4,
    type=int,
)
parser.add_argument(
    "-p",
    "--psl",
    help=".psl file from BLAT",
    required=True
)

parser.add_argument(
    "-i",
    "--cut_identity",
    help="Percentage Identity (fraction or percentage) for a hit to be counted",
    default=0.9,
    type=float
)

args = parser.parse_args()

cutoff = args.cutoff
bed = args.bed
psl=args.psl
cut_identity=args.cut_identity
if cut_identity<0:
    sys.exit(f"Cut Identity cannot be less that 0")
if cut_identity<=1:
    cut_identity*=100
if cut_identity>100:
    sys.exit(f"Cut Identity cannot be greater than 100%")

hit_count = defaultdict(int)

blat_headers=["matches", "misMatches", "repMatches", "nCount", 
      "qNumInsert", "qBaseInsert", "tNumInsert", "tBaseInsert", 
      "strand", "qName", "qSize", "qStart", "qEnd", "tName", 
      "tSize", "tStart", "tEnd", "blockCount", "blockSizes"]
with open(psl, "rt") as f:
    #try/except for empty input files
    try:
        # skip first 5 lines of psl file as it's headers
        for _ in range(5):
            next(f)
        for line in f:
            line = line.strip().split("\t")
            line_data=dict(zip(blat_headers,line))
            identity = float(line_data["matches"])/float(line_data["qSize"])*100
            if identity>= cut_identity:
                hit_count[line_data["qName"]] += 1
    except StopIteration:
        pass


with open(bed.replace(".bed", "_BLAT_filtered.bed"), "wt") as f_out:
    with open(bed, "rt") as f:
        for line in f:
            line = line.strip().split("\t")
            region = f"{line[0]}:{line[1]}-{line[2]}"
            if hit_count[region] <= cutoff:
                f_out.write("\t".join(line) + "\n")

