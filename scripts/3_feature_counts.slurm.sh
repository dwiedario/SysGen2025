#!/bin/bash
#SBATCH --job-name=featureCounts_mcortex
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=4G
#SBATCH --time=04:00:00
#SBATCH --output=log_star_align/featureCounts_%j.log

source ~/.bashrc
conda activate rnaseq

BASEDIR=/cluster/scratch/tkitak/brain_aging_rnaseq
STAR_DIR=/cluster/scratch/tkitak/star_align
GTF=/cluster/scratch/tkitak/ref/Mus_musculus.GRCm38.102.gtf

mkdir -p "${BASEDIR}/counts"

featureCounts \
  -T "${SLURM_CPUS_PER_TASK}" \
  -p -B -C \
  -a "${GTF}" \
  -o "${BASEDIR}/counts/mcortex_featureCounts.txt" \
  -t exon -g gene_id \
  ${STAR_DIR}/*.uniq.sorted.bam
