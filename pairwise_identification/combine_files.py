import argparse
import glob
import gzip
import os
import subprocess
import sys

import bx.align.axt
from natsort import natsorted
from tqdm import tqdm

parser = argparse.ArgumentParser()
parser.add_argument(
    "-d",
    "--directory",
    help="Directory or .net.axt files (or .net.axt.gz), split by chromosome to combine",
    required=True,
)
parser.add_argument(
    "-o",
    "--output",
    help="Name of output file, if not provided it will be based on the input file names",
)

if len(sys.argv) <= 1:
    sys.exit(parser.print_help())

args = parser.parse_args()

files = glob.glob(f"{args.directory}/*.net.axt*")

files = natsorted(files)
extension = set([file.split(".")[-1] for file in files])
if len(extension) > 1:
    sys.exit("Multiple file extensions found, please provide only one file extension")

extension = extension.pop()
if extension not in ("gz", "axt"):
    sys.exit(f"Invalid extension found: {extension}")

if extension == ".gz":
    open_func = gzip.open
else:
    open_func = open
outfile = args.output
if not outfile:
    outfile = os.path.join(os.path.dirname(files[0]), ".".join(files[0].split(".")[1:]))

if outfile in files:
    files.remove(outfile)

with open_func(outfile, "wt") as f_out:
    out_writer = bx.align.axt.Writer(f_out)
    for file in files:
        print(f"Processing: {file}")
        with open_func(file, "rt") as f:
            try:
                last_line = subprocess.check_output(["tail", "-4", file])
                last_record = int(last_line.decode().split("\n")[0].split(" ")[0])
            except UnicodeDecodeError:
                lines = f.read().splitlines()
                last_record = int(lines[-4].split("\n")[0].split(" ")[0])
                f.seek(0)
            axt_data = bx.align.axt.Reader(f)
            for record in tqdm(axt_data, total=last_record, leave=False):
                out_writer.write(record)
print(f"Output saved as: {outfile}")
