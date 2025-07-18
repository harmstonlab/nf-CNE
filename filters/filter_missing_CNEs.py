import glob
import sys
from collections import defaultdict
files = glob.glob("*combined_filtered.bed")

num_files = len(files)

if num_files == 0:
    print("No files found")
    sys.exit()

CE_seen_in_N_species=defaultdict(int)

for file in files:
    CE_in_species=set()
    with open(file,"rt") as f:
        for line in f:
            CNE_number=line.split("\t")[3]
            if CNE_number not in CE_in_species:
                CE_in_species.add(CNE_number)
                CE_seen_in_N_species[CNE_number]+=1

for file in files:
    with open(file) as f_in:
        with open(file.replace("combined_filtered.bed", "FINAL_filtered.bed"), "wt") as f_out:
            for line in f_in:
                CNE_number=line.split("\t")[3]
                if CE_seen_in_N_species[CNE_number] == num_files:
                    f_out.write("\t".join(line.strip().split("\t")[:5])+"\n")
