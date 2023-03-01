#!/bin/bash
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --export=ALL
#SBATCH --gres=gpu:v100-sxm2:1
#SBATCH --time=08:00:00
#SBATCH --mem=8G
#SBATCH --output=/work/donglab/ching.ki/ssRNA-MD/MD_sim/log/exec.%j.out

# usage: pipeline.sh <pdb> <conc> <ion>
# <pdb>  $1 path of pdb file
# <conc> $2 concentration of ion in mM
# <ion>  $3 type of ion to be added: Na+ or MG

SRC=/work/donglab/ching.ki/ssRNA-MD/MD_sim/src
DATA=/work/donglab/ching.ki/ssRNA-MD/VAE/data/gro
source $SRC/load_MD_tools.sh

# create unique directory to store intermediate files
#   ${1##*/} extract the name of the pdb file without the path to its directory
#   run_id can be used as a unique directory name since each sbtach runs simulation
#   for one pdb file for one ion in one concentration. Thus, each run have a unique
#   name so that concurrently running sbatch will not share the same directory
pdb_file=${1##*/}		# remove string leading last /
run_id=${pdb_file%.*}_$2_$3 	# remove string trailing last .
temp_dir=/scratch/ching.ki/$run_id
[ -d $temp_dir ] &&
	 rm -r $temp_dir &&
	 echo "removed old "$temp_dir" to create new one for this run"
mkdir "$temp_dir"
echo "Created new "$temp_dir
cd "$temp_dir"
 
# ./solvate.sh <pdb> <conc> <ion>
bash $SRC/solvate.sh $1 $2 $3 $SRC 
# python convert_prmtop-rst7_to_top-gro.py <amber-format-file-prefix>
python $SRC/convert_prmtop-rst7_to_top-gro.py rna 

# Run GROMACS
MDP_DIR=$SRC/md_params
# Energy minimization
# - assemble binary input using grompp with minim.mdp
# - carry out energy minimization
gmx grompp -f $MDP_DIR/minim.mdp -c rna.gro -p rna.top -o em.tpr
gmx mdrun -v -deffnm em 

# NVT ensemble equilibration
# - assemble binary input using grompp with nvt.mdp
# - carry out NVT ensemble equilibration
gmx grompp -f $MDP_DIR/nvt.mdp -c em.gro -r em.gro -p rna.top -o nvt.tpr
gmx mdrun -deffnm nvt

# NPT ensemble equilibration
# - assemble binary input using grompp with npt.mdp
# - carry out NPT ensemble equilibration
gmx grompp -f $MDP_DIR/npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt  -p rna.top -o npt.tpr
gmx mdrun -deffnm npt

# Production MD
# - assemble binary input using grompp with md.mdp
# - execute mdrun
gmx grompp -f $MDP_DIR/md.mdp -c npt.gro -t npt.cpt  -p rna.top -o md_0_1.tpr
gmx mdrun -deffnm md_0_1 -nb gpu

# Convert output to gro format for MDtraj python package to parse
# n_residues=$(python -c "import mdtraj as md; pdb=md.load('$1'); print(pdb.n_residues)")
gmx convert-trj -f md_0_1.xtc -s md_0_1.tpr -o $run_id".gro"

eval "$(conda shell.bash hook)"
conda activate py38
# Doesn't always work for some reason
#source /home/ching.ki/.bash_profile
#source activate py38
python -c "import mdtraj as md; pdb=md.load('$1'); t=md.load('$run_id.gro'); rna_idx=t.top.select(f'resi<{pdb.n_residues}'); rna=t.atom_slice(rna_idx); rna.save('$DATA/filtered_$run_id.gro')"

# Analyze trajectories for quality control
# TODO













# ARCHIVE
# (using mdtraj to process gro file instead)
# Remove unneeded atoms in gro output to save space
# grep -v "^[\ |0-9]*\(WAT\)\|\(Na\+\)\|\(Cl\-\)\|\(MG\)" $run_id".gro" > 
# 	$DATA/"filtered_"$run_id".gro"






