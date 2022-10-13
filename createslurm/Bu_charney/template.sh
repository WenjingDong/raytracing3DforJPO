#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16GB
#SBATCH --time=96:00:00
#SBATCH --output=slurm_%j.out

module purge
module load matlab/2020b
mfile=multraj3D
SRCDIR=$HOME/raytracing3D/source/charneyruns
INPDIR=$HOME/raytracing3D/main/Bu_charney/Bu16_a01/2f






