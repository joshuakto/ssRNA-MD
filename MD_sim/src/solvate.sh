#!/bin/bash
set -ex

# usage: solvate.sh
# ./solvate.sh <pdb> <conc> <ion>
# <pdb>  $1 name of odb file
# <conc> $2 concentration of ion in mM
# <ion>  $3 type of ion to be added: Na+ or MG
# <src>  $4 path to the source code used in this file
conc=$2
ion=$3
SRC=$4

## check what's the purpose pdb4amber -o $fn "$SCRIPTS/$1"

# Remove leap.log and volume.in if there exist from previous runs.
# Since these files will be used to infer the volume of the solution.
if [ -f leap.log ]; then
	rm leap.log
fi
if [ -f volume.in ]; then
	rm  volume.in
fi
# Unfortunately this is the only way I see how to get the volume
# of the simulation. Couldn't find how to load a prmtop and rst7
# file in tleap
cat > volume.in << EOF
source leaprc.RNA.OL3
source leaprc.water.opc
rna = loadpdb $1
addions rna Na+ 0
solvateoct rna OPCBOX 9.0
quit
EOF

tleap -f volume.in > /dev/null

ions=$(python $SRC/conc.py $conc < leap.log)
if [ $3 = "Na+" ]; then
	counter_ions=$ions
elif [ $3 = "MG" ]; then
	counter_ions=$((2*ions))
else
	echo "Unrecognized ion given: "$3
	exit 1
fi
rm leap.log volume.in

# Create the actual input file and redo everything
cat > tleap.in <<EOF
source leaprc.RNA.OL3
source leaprc.water.opc
rna = loadpdb $1
check rna
addions rna Na+ 0
solvateoct rna OPCBOX 9.0
addions rna $ion $ions
addions rna Cl- $counter_ions
saveamberparm rna rna.prmtop rna.rst7
quit
EOF

tleap -f tleap.in > tleap.out
