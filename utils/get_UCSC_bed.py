import argparse
import requests
import sys
import time
import io
import gzip
import re
import json

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-s", "--species", required=True, help="UCSC genome ID (e.g. mm10, hg38)")
    parser.add_argument("-t", "--track", required=True, help="UCSC track name (e.g. knownGene, rmsk)")
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
    return parser.parse_args()

def make_chrom_exclude(include_pattern, exclude_pattern):
    include_re = re.compile(include_pattern)
    exclude_re = re.compile(exclude_pattern)
    def _chrom_exclude(chrom):
        return not (include_re.match(chrom) and not exclude_re.search(chrom))
    return _chrom_exclude

def get_chromosomes(species,chrom_exclude):
    url = f"http://api.genome.ucsc.edu/list/chromosomes?genome={species}"
    r = requests.get(url)
    if r.status_code != 200:
        sys.exit(f"Failed to get chromosomes for {species}: {r.status_code}")
    return [chrom for chrom in r.json()["chromosomes"].keys() if not chrom_exclude(chrom)]

def detect_field(fields, candidates):
    for name in candidates:
        if name in fields:
            return name
    return None

# Some tracks (droPer1, ncbiRefSeq) have incorrect characters in them (some have a backspace?!?!)
# So if we can't make a JSON normally, sanitize the text and then make JSON
def sanitize_json_text(text):
    # Remove all ASCII control characters except tab and newline
    text = re.sub(r'[\x00-\x08\x0B-\x1F\x7F]', '', text)

    # Replace all backslashes not part of a valid JSON escape
    text = re.sub(r'\\(?![\"\\/bfnrtu])', r'\\\\', text)

    return text


def check_dataDownloadUrl(species, track):
    url = f"http://api.genome.ucsc.edu/getData/track?genome={species};track={track};maxItemsOutput=1;jsonOutputArrays=1"
    MAX_ATTEMPTS = 5
    for attempt in range(MAX_ATTEMPTS + 1):
        r = requests.get(url, timeout=(5, 5))
        if r.ok:
            r_json = r.json()
            if "dataDownloadUrl" in r_json.keys():
                track_data={"columnTypes":r_json["columnTypes"],track:[]}
                with requests.get(r_json["dataDownloadUrl"], stream=True) as data:
                    is_gzip = (r_json["dataDownloadUrl"].endswith(".gz") or data.headers.get("Content-Encoding") == "gzip")
                    byte_stream = gzip.GzipFile(fileobj=data.raw) if is_gzip else r.raw
                    text_data = io.TextIOWrapper(byte_stream, encoding="utf-8")
                    for line in text_data:
                        line=line.strip()
                        track_data[track].append(line.split("\t"))
                return track_data
        time.sleep((2 * attempt) * 0.5)
    return None

def fetch_track(species, track, chrom=None):
    if chrom:
        url = f"http://api.genome.ucsc.edu/getData/track?genome={species};track={track};chrom={chrom};maxItemsOutput=-1;jsonOutputArrays=1"
    else:
        url = f"http://api.genome.ucsc.edu/getData/track?genome={species};track={track};maxItemsOutput=-1;jsonOutputArrays=1"
    # Higher threshold if we are using chrom since we'll be making more requests so more likely to get rate limit
    MAX_ATTEMPTS = 5 if chrom else 2
    for attempt in range(MAX_ATTEMPTS + 1):
        print(f"{attempt=}")
        r = requests.get(url,timeout=60)
        if r.status_code == 200:
            try:
                return r.json()
            except (ValueError, json.JSONDecodeError):  
                try:
                    cleaned = sanitize_json_text(r.text)
                    return json.loads(cleaned)
                except Exception as e:
                    raise RuntimeError(f"JSON parse failed even after sanitization:\n{e}")
        if r.status_code == 206 and attempt == MAX_ATTEMPTS:
            try:
                return r.json() if not chrom else None
            except (ValueError, json.JSONDecodeError):
                cleaned = sanitize_json_text(r.text)
                return json.loads(cleaned)
        time.sleep((2 ** attempt) * 0.5)
    return None


def write_track_bed(data, track, bed_handle, chrom_exclude):
    column_map = {col["name"]: i for i, col in enumerate(data["columnTypes"][track])}

    chrom_field = detect_field(column_map, ["chrom", "genoName"])
    start_field = detect_field(column_map, ["chromStart", "genoStart","cdsStart", "txStart","thickStart"])
    end_field = detect_field(column_map, ["chromEnd", "genoEnd","cdsEnd", "txEnd","thickEnd"])

    cds_start = detect_field(column_map, ["cdsStart", "thickStart"])
    if cds_start:
        cds_start = column_map.get(cds_start)
    cds_end = detect_field(column_map, ["cdsEnd", "thickEnd"])
    if cds_end:
        cds_end = column_map.get(cds_end)
    exon_starts = column_map.get("exonStarts")
    exon_ends = column_map.get("exonEnds")
    exon_frames = column_map.get("exonFrames")

    chrom_starts = column_map.get("chromStarts")
    block_sizes = column_map.get("blockSizes")
    if chrom_field is None or start_field is None or end_field is None:
        raise RuntimeError("Missing required fields for BED output.")

    chrom_idx = column_map[chrom_field]
    start_idx = column_map[start_field]
    end_idx = column_map[end_field]

    records = data[track].values() if isinstance(data[track], dict) else data[track]
    for group in records:
        if not group:
            continue
        group = group if isinstance(group[0], list) else [group]
        for rec in group:
            chrom = rec[chrom_idx]
            if chrom_exclude(chrom):
                continue
            if exon_frames is not None and exon_starts is not None and exon_ends is not None:
                starts = list(map(int, rec[exon_starts].rstrip(",").split(",")))
                ends = list(map(int, rec[exon_ends].rstrip(",").split(",")))
                frames = list(map(int, rec[exon_frames].rstrip(",").split(",")))
                for start, end, frame in zip(starts, ends, frames):
                    if frame != -1:
                        if int(end)>int(rec[cds_end]):
                            end=int(rec[cds_end])
                        if int(start)<int(rec[cds_start]):
                            start=int(rec[cds_start])
                        bed_handle.write(f"{chrom}\t{start}\t{end}\n")

            elif chrom_starts is not None and block_sizes is not None  and exon_frames is not None and cds_start is not None and cds_end is not None:
                chrom_start_val = rec[start_idx]
                starts = list(map(int, rec[chrom_starts].rstrip(",").split(",")))
                sizes = list(map(int, rec[block_sizes].rstrip(",").split(",")))
                frames = list(map(int, rec[exon_frames].rstrip(",").split(",")))
                for offset, size, frame in zip(starts, sizes,frames):
                    if frame != -1:
                        start = int(chrom_start_val) + int(offset)
                        end = int(start) + int(size)
                        if int(end)>int(rec[cds_end]):
                            end=int(rec[cds_end])
                        if int(start)<int(rec[cds_start]):
                            start=int(rec[cds_start])
                        bed_handle.write(f"{chrom}\t{start}\t{end}\n")

            elif cds_start is not None and cds_end is not None and exon_starts is not None and exon_ends is not None:
                cds_start_val = rec[cds_start]
                cds_end_val = rec[cds_end]
                starts = list(map(int, rec[exon_starts].rstrip(",").split(",")))
                ends = list(map(int, rec[exon_ends].rstrip(",").split(",")))
                for start, end in zip(starts, ends):
                    coding_start = max(int(start), int(cds_start_val))
                    coding_end = min(int(end), int(cds_end_val))
                    if coding_start < coding_end:
                        bed_handle.write(f"{chrom}\t{coding_start}\t{coding_end}\n")


            else:
                start = rec[start_idx]
                end = rec[end_idx]
                bed_handle.write(f"{chrom}\t{start}\t{end}\n")

def main():
    args = parse_args()
    chrom_exclude = make_chrom_exclude(args.chr_include, args.chr_exclude)
    chromosomes = get_chromosomes(args.species, chrom_exclude)
    outname = f"{args.species}_{args.track}.bed"

    with open(outname, "wt") as bed:
        #First try and download entire thing directly with dataDownloadUrl
        data = check_dataDownloadUrl(args.species, args.track)
        #If that doesn't work retrieve records from API through API
        if not data:
            data = fetch_track(args.species, args.track)
        if data is None:
            raise RuntimeError(f"UCSC fetch failed or returned truncated data for {args.species} - {args.track}.")
        if "maxItemsLimit" not in data.keys():
            print(f"Downloading entire {args.track}")
            write_track_bed(data, args.track, bed, chrom_exclude)
        # If not (we hit max items of 1M per request) fall back to getting per chr
        # Species like droPer1 have many chroms/scaffolds but small size. So with rate limiting would take hours
        # If we did the entire thing per chr, or seconds if combined. But in other cases (e.g. hg38 rmsk)
        # The entire thing is too big so we MUST do it per chrom
        else:
            for chrom in chromosomes:
                print(f"Downloading {args.track} from {chrom}")
                data = fetch_track(args.species, args.track, chrom)
                if data is None:
                    raise RuntimeError(f"UCSC fetch failed or returned truncated data for {chrom}.")
                write_track_bed(data, args.track, bed, chrom_exclude)

    print(f"\nComplete. BED file written: {outname}")

if __name__ == "__main__":
    main()
