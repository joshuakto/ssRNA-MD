#!/usr/bin/env python
import re
import math
import sys
import argparse

parser = argparse.ArgumentParser(description="Calculate the number of Na+ and Cl- atoms to add from a leap.log as stdin")
parser.add_argument("conc", type=float, help="Target concentration in mM")
args = parser.parse_args()

matches = re.findall("Volume: (\d+\.\d+)", sys.stdin.read())
if len(matches) != 1:
    print(f"Expected exactly 1 Volume entry in leap.log found {len(matches)}")
    sys.exit(1)

volume = float(matches[0])
L = volume * 1E-27 # since tleap gives volume in cubic Amstrong, which is 1e-27 liter
conc = args.conc
atoms_per_L = (conc * 6.022E23) / 1E3
print(math.ceil(L * atoms_per_L))
