#!/bin/bash
#SBATCH --job-name=star_index_mm10
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=6G
#SBATCH --time=08:00:00
#SBATCH --output=log/star_index_%j.log

source ~/.bashrc
conda activate rnaseq

GENOME_FASTA=/cluster/scratch/tkitak/ref/Mus_musculus.GRCm38.dna.primary_assembly.fa
GTF_FILE=/cluster/scratch/tkitak/ref/Mus_musculus.GRCm38.102.gtf
INDEX_DIR=/cluster/scratch/tkitak/star_index

mkdir -p "${INDEX_DIR}"

STAR \
  --runThreadN "${SLURM_CPUS_PER_TASK}" \
  --runMode genomeGenerate \
  --genomeDir "${INDEX_DIR}" \
  --genomeFastaFiles "${GENOME_FASTA}" \
  --sjdbGTFfile "${GTF_FILE}" \
  --sjdbOverhang 99

echo "STAR index build complete"
