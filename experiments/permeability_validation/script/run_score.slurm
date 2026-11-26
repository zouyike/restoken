#!/bin/bash
#SBATCH --job-name=pm_2
#SBATCH --partition=gpu_cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH --time=00:30:00
#SBATCH --output=/scratch/genesis/NCAA_tokenization/experiments/permeability_validation/output/score_%j.log
#SBATCH --error=/scratch/genesis/NCAA_tokenization/experiments/permeability_validation/output/score_%j.log

export HF_HUB_OFFLINE=1
export TRANSFORMERS_OFFLINE=1
export TOKENIZERS_PARALLELISM=false
export PYTHONUNBUFFERED=1

PY=/public/home/genesis/miniconda3/envs/esm3/bin/python
SCRIPT=/scratch/genesis/NCAA_tokenization/experiments/permeability_validation/script/score_generations.py

echo "host=$(hostname) start=$(date)"
$PY $SCRIPT
echo "done=$(date)"
