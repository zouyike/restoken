#!/bin/bash
#SBATCH --job-name=cv_sp
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=/scratch/genesis/NCAA_tokenization/experiments/dd_revision/perm_splits/prep_%j.log
#SBATCH --error=/scratch/genesis/NCAA_tokenization/experiments/dd_revision/perm_splits/prep_%j.log

export PYTHONUNBUFFERED=1
/public/home/genesis/miniconda3/envs/ai/bin/python /scratch/genesis/NCAA_tokenization/experiments/dd_revision/perm_splits/prepare_splits_and_baselines.py
