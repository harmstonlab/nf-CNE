import argparse
import sys

import requests

parser = argparse.ArgumentParser()
parser.add_argument(
    "-s",
    "--species",
    help="Species name to query against UCSC api",
    required=True,
)

if len(sys.argv) <= 1:
    sys.exit(parser.print_help())

args = parser.parse_args()
species = args.species
search_url = f"http://api.genome.ucsc.edu/list/tracks?genome={species};trackLeavesOnly=1"
response=requests.get(search_url)
if response.status_code!=200:
    try:
        message=response.json()["error"]
    except:
        message="NOT PROVIDED"
    sys.exit(f"Error Accessing UCSC API. Status code: {response.status_code}. Error cause: {message}")

else:
    json=response.json()
    species_data=None
    try:
        species_data=json[species]
    except KeyError:
        for key in json.keys():
            if type(json[key])==dict:
                species_data=json[key]
                break
    if species_data is None:
        sys.exit(f"Error: No data found for species '{species}' in UCSC API response.")
    for track in species_data:
        if "compositeContainer" in species_data[track].keys():
            continue
        elif "bed" in species_data[track]["type"].lower():
            print(track.strip())
        elif "genepred" in species_data[track]["type"].lower():
            print(track.strip())
        elif "rmsk" in species_data[track]["type"].lower():
            print(track.strip())