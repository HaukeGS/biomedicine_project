#!/bin/bash
#SBATCH --job-name=scFoundation
#SBATCH --partition=gpu

#SBATCH --output=output.log
#SBATCH --error=fehler.log
#SBATCH --time=02:00:00
#SBATCH --gres=gpu:a100m40:1
#SBATCH --mem=256G
#SBATCH --cpus-per-task=12

module load Miniforge3
conda activate scFoundation

python script.py