#!/bin/bash
#SBATCH --job-name=pm_1
#SBATCH --partition=gpu_cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH --time=08:00:00
#SBATCH --output=/scratch/genesis/NCAA_tokenization/experiments/permeability_validation/output/train_%j.log
#SBATCH --error=/scratch/genesis/NCAA_tokenization/experiments/permeability_validation/output/train_%j.log

export HF_HUB_OFFLINE=1
export TRANSFORMERS_OFFLINE=1
export TOKENIZERS_PARALLELISM=false
export PYTHONUNBUFFERED=1

PY=/public/home/genesis/miniconda3/envs/esm3/bin/python
SCRIPT=/scratch/genesis/NCAA_tokenization/experiments/permeability_validation/script/train_perm.py

echo "host=$(hostname) gpu=$CUDA_VISIBLE_DEVICES start=$(date)"
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader
$PY $SCRIPT --epochs 40 --bs 32 --lr 2e-5 --maxlen 512 --patience 6
echo "done=$(date)"
