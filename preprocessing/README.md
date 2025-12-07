# systems-genomics-tutorial

# How I did Download
prefetch --option-file SRR_Acc_List.txt -O raw_data/

cd raw_data
for sra in SRR*/*.sra; do
    fasterq-dump --threads 8 --progress "$sra" --outdir ../fastq/
done

# Reference 

wget -c https://ftp.ensembl.org/pub/release-113/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
wget -c https://ftp.ensembl.org/pub/release-113/gtf/mus_musculus/Mus_musculus.GRCm39.113.gtf.gz


# Generated raw count x sample count matrix from RSEM output:
rsem-generate-data-matrix *.genes.results > gene_counts_matrix.txt