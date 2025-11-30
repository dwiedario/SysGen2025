#!/bin/bash
#SBATCH --job-name=trim_mcortex
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G
#SBATCH --time=04:00:00
#SBATCH --array=0-21
#SBATCH --output=/cluster/scratch/tkitak/log/trim_%A_%a.log

source ~/.bashrc
conda activate rnaseq

FASTQ_DIR=/cluster/scratch/tkitak/brain_aging_rnaseq/fastq
TRIM_DIR=/cluster/scratch/tkitak/brain_aging_rnaseq/trimmed

mkdir -p "${TRIM_DIR}"

SAMPLES_FILE=/cluster/scratch/tkitak/systems-genomics-tutorial/SRR_Acc_List.txt
SAMPLE=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "${SAMPLES_FILE}")

R1=${FASTQ_DIR}/${SAMPLE}_1.fastq
R2=${FASTQ_DIR}/${SAMPLE}_2.fastq

cd "${TRIM_DIR}"

trim_galore \
  --paired \
  --length 20 \
  --phred33 \
  -q 30 \
  "${R1}" "${R2}"
