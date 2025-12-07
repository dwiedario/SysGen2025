#!/bin/bash
#SBATCH --job-name=fastqc           # Job name
#SBATCH --ntasks=1                  # One task
#SBATCH --cpus-per-task=4           # Number of CPU cores
#SBATCH --mem-per-cpu=2G            # Memory per core
#SBATCH --time=04:00:00             # Walltime (hh:mm:ss)
#SBATCH --output=logs/fastqc_%j.log # Log file

WORKDIR=/cluster/scratch/tkitak/brain_aging_rnaseq
cd "${WORKDIR}"

mkdir -p logs
mkdir -p fastqc_reports

source ~/.bashrc
conda activate rnaseq

fastqc -t "${SLURM_CPUS_PER_TASK}" fastq/*.fastq --outdir fastqc_reports
