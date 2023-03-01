#!/bin/bash
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --export=ALL
#SBATCH --gres=gpu:1
#SBATCH --time=00:20:00
#SBATCH --mem=2G
#SBATCH --output=/work/donglab/ching.ki/ssRNA-MD/VAE/log/process_traj/exec.%j.out

# usage: process_trajectory.sh <filename> <gro_dir> <out_dir>
# <filename>  $1 name of the trajectory to be processed, i.e., filtered_sasdfb9_<model>_<conc>_<ion> 
# <gro_dir>   $2 path of the directory storing gro input file storing the trajectory
# <out_dir>   $3 path to directory to store the output numpy compressed(npz) file containing the relative distances between atoms

eval "$(conda shell.bash hook)"
conda activate py38

python process_trajectory.py -f $1 -gro_dir $2 -out_dir $3
