#!/bin/bash 
#SBATCH --partition=shared
#SBATCH -A uoo104
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=100G
#SBATCH --time 24:00:00
#SBATCH --job-name Pmod_proc
#SBATCH --output Pmod_proc-%J.log
#SBATCH --mail-user=jwelshka@uoregon.edu
#SBATCH --mail-type=begin  # email me when the job starts
#SBATCH --mail-type=end    # email me when the job finishes

# this script will turn a GROMACS .g96 file into numpy array and save it as a .npy file

BD=$PWD
cd $BD

# loading GROMACS may have slightly different prerequisites on different platforms

# 1.

module load intel
module load cpu/0.15.4
module load gcc/10.2.0
module load openmpi/4.0.4
module load gromacs/2020.4
module list

# 2.

echo "3" | gmx_mpi trjconv -f after_rot_example.xtc -s SP_dsAXA.pdb -n SP_index.ndx -o P_traj.g96

# 3.
gmx_mpi sasa -f after_rot_example.xtc -s SP_dsAXA.pdb -n SP_index.ndx -or P_resarea.xvg -dt 1000 -surface 1 -output 2 3

# 4.

sed '/BOX/, +1 d' P_traj.g96 | sed '/TITLE/, +1 d' | awk 'NF==3' > Ptmp

echo "Phosphate trajectory loaded into txt file!"

python -c 'import numpy as np; np.save("P_traj.npy",np.loadtxt("Ptmp"))'

echo "Phosphate trajectory converted to .npy file!"

rm -rf Ptmp

# 5.

mkdir $BD/tcfrun

python process_Pmodel.py

# 6.

cp tcf.sh $BD/tcfrun/

cp tcfint_DUMMY.f $BD/tcfrun/

cp tcf_DUMMY.pbs $BD/tcfrun/

cp move.sh $BD/tcfrun/

cp protname.txt $BD/tcfrun

cd $BD/tcfrun/

sh $BD/tcfrun/tcf.sh

cd $BD

# 7.

touch nmol.dat

echo "2" > nmol.dat

cp mode_analysis/barriers_kcal.dat $BD

gfortran m1calc.f -o m1calc.exe

./m1calc.exe

# 8.

mkdir numpy_arrays

mv *.npy numpy_arrays

mkdir $BD/tcfdata

mv m1_* $BD/tcfdata

cd $BD

# 9.

touch LE4DNA_OUTPUT

echo "If there are no errors in the .log file then the LE4DNA analysis is complete." >> LE4DNA_OUTPUT

echo "You may now run move.sh in the tcfrun directory, provided the tcfs have been calculated. The tcfs of bond autocorrelation may take several hours to calculate." >> LE4DNA_OUTPUT




