import mdtraj as md
import numpy as np
np.set_printoptions(precision=4)
import pandas as pd
import math
from tqdm import tqdm
import concurrent.futures
import argparse


def parse_args():
    parser = argparse.ArgumentParser(
        description = 'Process given trajectory and return an array containing the relative distances between atoms for each frame in the trajectory.')
    parser.add_argument('-f', '--filename', type=str,
                        help='name of the trajectory to be processed, i.e., filtered_sasdfb9_<model>_<conc>_<ion>')
    parser.add_argument('-gro_dir', type=str, default="../data/gro",
                        help='path of the directory storing gro input file storing the trajectory.')
    parser.add_argument('-out_dir', type=str, default="../data/distance_npz",
                        help='path to directory to store the output numpy compressed(npz) file containing the relative distances between atoms.')
    args = parser.parse_args()
    return args


def compute_and_pack_distance_in_array_frames(rna, pair_grid):
    # in distances_frames, the axes were (number of columns, number of frames, number of rows)
    # where number of columns and rows are both equal to the number of atoms
    distances_by_columns = np.array([md.compute_distances(rna, column) for column in pair_grid])
    # in axes_ordered_distances_frames, the axes are 
    # (number of frames, number of columns, number of rows)
    axes_ordered_distances_frames = distances_by_columns.transpose(1,0,2)
    packed_distances_frames = np.apply_along_axis(
        lambda arr: np.array_split(arr, len(arr)), 2, axes_ordered_distances_frames)
    return(packed_distances_frames)

# Given an 3D array (n_atom x n_atom x 1), n_atom is the number of atom in one batch of residue and 
# may vary depending on the number of atoms in each residue; return a 100x100x1 array either by 
# truncating or padding the given array
def trunc_or_pad_distance_array(distance_array, maxlen=100):
    if distance_array.shape[0] != distance_array.shape[1]:
        raise Exception("the first and second dimension of the given distance array does not match")
    n_atom = distance_array.shape[0]
    if n_atom < maxlen:
        pad_n_pre = int((maxlen-n_atom)/2)
        pad_n_post = (maxlen-n_atom)/2
        if pad_n_pre != pad_n_post:
            pad_n_post = math.ceil(pad_n_post)
        distance_array = np.pad(
            distance_array, 
            ((pad_n_pre, pad_n_post), 
             (pad_n_pre, pad_n_post),
             (0, 0)))
    elif n_atom > maxlen:
        trunc_n_pre = int((n_atom-maxlen)/2)
        trunc_n_post = (n_atom-maxlen)/2
        if trunc_n_pre != trunc_n_post:
            trunc_n_post = math.ceil(trunc_n_post)
        end_idx = int(n_atom - trunc_n_post)
        # truncate first dimension
        distance_array = distance_array[trunc_n_pre:end_idx]
        # truncate second dimension
        distance_array = np.array([col[trunc_n_pre:end_idx] for col in distance_array])    
    return distance_array

def preprocess_atom_set(rna, atoms_idx):
    print(f"Processing {atoms_idx}.")
    pair_meshgrid = np.meshgrid(atoms_idx, np.transpose(atoms_idx))
    pair_grid = np.dstack((pair_meshgrid[0], pair_meshgrid[1]))
    # TODO can be optimized to only calculate distance once for each pair but the computation quick so I am skipping for now
    pair_distance_array_frames = compute_and_pack_distance_in_array_frames(rna, pair_grid)
    formatted_pair_distance_array_frames = np.array(
        [trunc_or_pad_distance_array(frame, maxlen = 100) 
         for frame in pair_distance_array_frames])
    print(f"Processed {atoms_idx}.")
    return formatted_pair_distance_array_frames

def main(filename, gro_dir, out_dir):
    gro_file_path = f"{gro_dir}/{filename}.gro"
    outfile_path = f"{out_dir}/{filename}.npz"
    print(f"Reading trajectory from: {gro_file_path}")
    rna = md.load(gro_file_path)
    print(f"Trajectory loaded from: {gro_file_path}")
    n_residues = rna.top.n_residues
    n_atoms = rna.top.n_atoms
    
    
    ## First trial: select atoms by groups of residue, every 3 consecutive residues as a group
    ## Issue: Does not show great clustering by meaningful variables, 
    ## maybe because it only captures local information, 
    ## which might not depend on variables such as type of ion and ionc concentration. 
    # batch_residue_idx = [(i, i+3) for i in range(0,n_residues,3)]
    # atoms_idx_list = [rna.top.select(f'resi >= {batch_start} and resi < {batch_end}')
    #                 for batch_start, batch_end in batch_residue_idx]
    
    ## Second trial, select 100 atoms in the model distributed across the entire ssRNA molecule
    interval = math.floor(n_atoms/100)
    atom_ids = [i for i in range(0, n_atoms, interval)][0:100]
    atoms_idx_list = [atom_ids]
    
    n_atoms_idx_set = len(atoms_idx_list)
    print(f"{n_residues} residues found in trajectory")
    # print(f"Processing atoms from residues in the following batches: {batch_residue_idx}")
    print(f"Extracting atoms with idx: {atom_ids}")
    # # For some reason concurrent.futures does not work here
    # with concurrent.futures.ProcessPoolExecutor() as executor:
    #     dataset_array_list = executor.map(
    #         lambda aset: preprocess_atom_set(rna, aset), 
    #         atoms_idx_list)
    dataset_array_list = [
        preprocess_atom_set(rna, atoms_idx_list[i])
        for i in tqdm(range(n_atoms_idx_set))]
    dataset_array = np.concatenate(dataset_array_list)
    print("Processing for all atom sets done.")
    np.savez_compressed(outfile_path, dataset_array)
    print(f"dataset_array saved to {outfile_path}")
    return dataset_array

if __name__ == "__main__":
    args = parse_args()
    main(args.filename, args.gro_dir, args.out_dir)