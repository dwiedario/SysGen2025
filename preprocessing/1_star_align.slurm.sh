#!/bin/bash
#SBATCH --job-name=star_align_mcortex
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=5G
#SBATCH --time=04:00:00
#SBATCH --array=0-21
#SBATCH --output=log_star_align/star_%A_%a.log

source ~/.bashrc
conda activate rnaseq

TRIM_DIR=/cluster/scratch/tkitak/brain_aging_rnaseq/trimmed
STAR_DIR=/cluster/scratch/tkitak/star_align
GENOME_DIR=/cluster/scratch/tkitak/star_index

mkdir -p "${STAR_DIR}"

SAMPLES_FILE=/cluster/scratch/tkitak/systems-genomics-tutorial/SRR_Acc_List.txt
SAMPLE=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "${SAMPLES_FILE}")

R1=${TRIM_DIR}/${SAMPLE}_1_val_1.fq
R2=${TRIM_DIR}/${SAMPLE}_2_val_2.fq

OUT_PREFIX=${STAR_DIR}/${SAMPLE}.

STAR \
  --runThreadN "${SLURM_CPUS_PER_TASK}" \
  --genomeDir "${GENOME_DIR}" \
  --readFilesIn "${R1}" "${R2}" \
  --outFileNamePrefix "${OUT_PREFIX}" \
  --outFilterMultimapNmax 1 \
  --quantMode TranscriptomeSAM \
  --outSAMtype BAM SortedByCoordinate

echo "Alignment done"
