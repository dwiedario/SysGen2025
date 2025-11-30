# systems-genomics-tutorial

# How I did Download
prefetch --option-file SRR_Acc_List.txt -O raw_data/

cd raw_data
for sra in SRR*/*.sra; do
    fasterq-dump --threads 8 --progress "$sra" --outdir ../fastq/
done
