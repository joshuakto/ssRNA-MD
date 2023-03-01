# use conda environment py381
import parmed as pmd
import argparse

parser = argparse.ArgumentParser(description='Convert amber files into gromacs format to carry out MD simulation in GROMACS.')
parser.add_argument('mol', type=str, nargs=1,
                    help='molecule for which the files are converted')
args = parser.parse_args()
parm = pmd.load_file(f"{args.mol[0]}.prmtop", f"{args.mol[0]}.rst7")
parm.save(f"{args.mol[0]}.top", format="gromacs")
parm.save(f"{args.mol[0]}.gro")
