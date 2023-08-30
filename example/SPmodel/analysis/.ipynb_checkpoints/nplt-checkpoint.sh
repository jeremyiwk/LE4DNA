#!/bin/bash 
#SBATCH --partition=shared
#SBATCH -A uoo104
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=100G
#SBATCH --time 24:00:00
#SBATCH --job-name loadtxt
#SBATCH --output loadtxt-%J.log
#SBATCH --mail-user=jwelshka@uoregon.edu
#SBATCH --mail-type=begin  # email me when the job starts
#SBATCH --mail-type=end    # email me when the job finishes

# this script will turn a GROMACS .g96 file into numpy array and save it as a .npy file

BD=$PWD
cd $BD

sed '/BOX/, +1 d' S_traj.g96 | sed '/TITLE/, +1 d' | awk 'NF==3' > Stmp

echo "1:done!"

sed '/BOX/, +1 d' P_traj.g96 | sed '/TITLE/, +1 d' | awk 'NF==3' > Ptmp

echo "2:done!"

python -c 'import numpy as np; np.save("P_traj.npy",np.loadtxt("Ptmp"))'

echo "3:done!"

python -c 'import numpy as np; np.save("S_traj.npy",np.loadtxt("Stmp"))'

echo "4:done!"
