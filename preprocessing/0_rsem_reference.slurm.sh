#!/bin/bash
#SBATCH --job-name=rsem_prep_GRCm39
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=6G
#SBATCH --time=04:00:00
#SBATCH --output=log/rsem_prep_%j.log

source ~/.bashrc
conda activate rnaseq

GENOME_FASTA=/cluster/scratch/tkitak/ref/Mus_musculus.GRCm39.dna.primary_assembly.fa
GTF_FILE=/cluster/scratch/tkitak/ref/Mus_musculus.GRCm39.113.gtf
RSEM_REF_BASE=/cluster/scratch/tkitak/rsem_ref/Mus_musculus.GRCm39.113

mkdir -p "$(dirname "${RSEM_REF_BASE}")"

rsem-prepare-reference \
  --gtf "${GTF_FILE}" \
  -p "${SLURM_CPUS_PER_TASK}" \
  "${GENOME_FASTA}" \
  "${RSEM_REF_BASE}"

echo "RSEM reference built at ${RSEM_REF_BASE}"
