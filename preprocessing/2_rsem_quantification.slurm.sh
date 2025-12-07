#!/bin/bash
#SBATCH --job-name=rsem_quant_mcortex
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=5G
#SBATCH --time=04:00:00
#SBATCH --array=0-21
#SBATCH --output=log_rsem/rsem_%A_%a.log

source ~/.bashrc
conda activate rnaseq

STAR_DIR=/cluster/scratch/tkitak/star_align
RSEM_REF_BASE=/cluster/scratch/tkitak/rsem_ref/Mus_musculus.GRCm39.113
RSEM_OUT_DIR=/cluster/scratch/tkitak/rsem_quant
SAMPLES_FILE=/cluster/scratch/tkitak/systems-genomics-tutorial/SRR_Acc_List.txt

mkdir -p "${RSEM_OUT_DIR}"

SAMPLE=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "${SAMPLES_FILE}")

IN_BAM=${STAR_DIR}/${SAMPLE}.Aligned.toTranscriptome.out.bam
OUT_PREFIX=${RSEM_OUT_DIR}/${SAMPLE}

rsem-calculate-expression \
  --bam \
  --no-bam-output \
  --paired-end \
  -p "${SLURM_CPUS_PER_TASK}" \
  --forward-prob 0.5 \
  "${IN_BAM}" \
  "${RSEM_REF_BASE}" \
  "${OUT_PREFIX}"

echo "RSEM quantification done for ${SAMPLE}"
