#!/bin/bash
#SBATCH --job-name=filter_sort_mcortex
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G
#SBATCH --time=04:00:00
#SBATCH --array=0-21
#SBATCH --output=log_star_align/filter_sort_%A_%a.log

source ~/.bashrc
conda activate rnaseq

STAR_DIR=/cluster/scratch/tkitak/star_align
SAMPLES_FILE=/cluster/scratch/tkitak/systems-genomics-tutorial/SRR_Acc_List.txt

SAMPLE=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "${SAMPLES_FILE}")

SAM=${STAR_DIR}/${SAMPLE}.Aligned.out.sam
BAM_UNIQ_SORTED=${STAR_DIR}/${SAMPLE}.uniq.sorted.bam

# Keep only header lines and reads with NH:i:1 (uniquely mapped), then convert to BAM and sort
samtools view -h "$SAM" \
  | awk '$1 ~ /^@/ || $0 ~ /NH:i:1/' \
  | samtools view -b - \
  | samtools sort -@ "${SLURM_CPUS_PER_TASK}" -o "$BAM_UNIQ_SORTED"

samtools index "$BAM_UNIQ_SORTED"

echo "Done filter sort BAM"