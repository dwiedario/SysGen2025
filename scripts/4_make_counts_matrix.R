#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(rtracklayer)
})

fc_file  <- "/cluster/scratch/tkitak/brain_aging_rnaseq/counts/mcortex_featureCounts.txt"
gtf_file <- "/cluster/scratch/tkitak/ref/Mus_musculus.GRCm38.102.gtf"
out_file <- "/cluster/scratch/tkitak/brain_aging_rnaseq/counts/mcortex_counts_protein_coding.tsv"

## 1. Read featureCounts output
fc <- read_tsv(
  fc_file,
  comment = "#",
  col_types = cols()
)

# annotation columns + counts
counts <- as.matrix(fc[, 7:ncol(fc)])
rownames(counts) <- fc$Geneid

## 2. Clean column names to SRR IDs
old_names   <- colnames(counts)
clean_names <- sub("\\.uniq\\.sorted\\.bam$", "", basename(old_names))
colnames(counts) <- clean_names

## 3. Load GTF, get gene-level annotation
gtf <- import(gtf_file)
gtf_genes <- gtf[gtf$type == "gene"]
m <- mcols(gtf_genes)

biotype_col <- if ("gene_biotype" %in% colnames(m)) "gene_biotype" else "gene_type"

annot <- as.data.frame(m[, c("gene_id", "gene_name", biotype_col)])
colnames(annot)[3] <- "gene_biotype"

# Align to count matrix rows
annot <- annot[match(rownames(counts), annot$gene_id), ]

## 4. Keep only protein-coding genes
is_pc <- !is.na(annot$gene_biotype) & annot$gene_biotype == "protein_coding"

counts_pc <- counts[is_pc, , drop = FALSE]
gene_symbols <- annot$gene_name[is_pc]

# collapse duplicated symbols by summing
if (any(duplicated(gene_symbols))) {
  counts_pc <- rowsum(counts_pc, group = gene_symbols)
} else {
  rownames(counts_pc) <- gene_symbols
}

cat("Rows after protein-coding filter:", nrow(counts_pc), "\n")

## 5. Write simple matrix (gene symbol + SRR columns)
out_df <- data.frame(
  gene = rownames(counts_pc),
  counts_pc,
  check.names = FALSE
)

write_tsv(out_df, out_file)
