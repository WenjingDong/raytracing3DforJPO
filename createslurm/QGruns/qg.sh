#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=8GB
#SBATCH --time=124:00:00
#SBATCH --mail-type=END
#SBATCH --mail-user=wd583@nyu.edu
#SBATCH --output=slurm_%j.out

module purge
module load matlab/2020b
mfile=multraj3D_scaled
SRCDIR=$HOME/raytracing3D/source/QGruns
INPDIR=$HOME/raytracing3D/main/QGN100/scaled/2f






