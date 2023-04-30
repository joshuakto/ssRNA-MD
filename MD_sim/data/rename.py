import json
import argparse
import re
import sys
import os
from pprint import pprint

parser = argparse.ArgumentParser(description="Rename RNA atoms for rU30 and rA30")
parser.add_argument("pdb", type=str, help="The pdb file to rename") 

args = parser.parse_args()

def panic(string):
    print(string, file=sys.stderr)
    sys.exit(1)

fn_re = re.compile("(\w+)_fit(\d+)_model(\d+)\.pdb")
head, tail = os.path.split(args.pdb)
match = fn_re.match(tail)
if not match:
    panic("File must match the syntax '<code>_fit<n>_model<m>.pdb'")

code = match.group(1)
fit = match.group(2)
model = match.group(3)

with open("mapping.json", "r") as f:
    config = json.load(f)

if code not in config:
    panic(f"Code not recognized: {code}")

# 1 Indexed list of lines to change
def get_lines(which):
    if which == "A":
        first_size = 19
    elif which == "U":
        first_size = 17
    else:
        panic("Bad get lines")
    
    offset = first_size + 1
    res_size = first_size + 3

    yield from range(1, 10)              # First sugar atoms
    for i in range(30):
        yield res_size*i   + 1  + offset # Oxygen
        start = res_size*i + 4  + offset # Sugar start
        end   = res_size*i + 11 + offset # Sugar end
        yield from range(start, end + 1)

outfn = open(f"{code.lower()}_m{model}.pdb", "w")
acid = config[code]["acid"]
update = set(get_lines(acid))
with open(args.pdb, "r") as f:
    for i, line in enumerate(f):
        line = line.replace(f"r{acid}", f" {acid}")
        if i + 1 in update:
            name = line[12:16]
            new = name
            if not name.endswith("'"):
                new = name.strip() + "'"
                new = new.rjust(4, " ")
                line = line.replace(name, new)
        print(line, file=outfn, end="")

outfn.close()
