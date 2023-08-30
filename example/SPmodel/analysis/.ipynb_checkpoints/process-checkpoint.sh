#!/bin/bash 
#SBATCH --partition=shared
#SBATCH -A uoo104
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=200G
#SBATCH --time 24:00:00
#SBATCH --job-name process-%j
#SBATCH --output process-%j.log
#SBATCH --mail-user=jwelshka@uoregon.edu
#SBATCH --mail-type=begin  # email me when the job starts
#SBATCH --mail-type=end    # email me when the job finishes


BD=$PWD
cd $BD

python process.py
