library(rtracklayer)

## 1. Load count matrix from RSEM
counts <- read.table(
  "gene_counts_matrix.txt",
  header = TRUE,
  row.names = 1,
  check.names = FALSE,
  quote = "\""
)

# Clean column names
colnames(counts) <- sub("\\.genes\\.results$", "", colnames(counts))

## 2. Load annotation from GTF
gtf_file <- "/cluster/scratch/tkitak/ref/Mus_musculus.GRCm39.113.gtf"
gtf <- import(gtf_file)
genes <- gtf[gtf$type == "gene"]

anno <- data.frame(
  gene_id   = genes$gene_id,
  biotype   = as.character(genes$gene_biotype),
  gene_name = as.character(genes$gene_name),
  stringsAsFactors = FALSE
)

# Align annotation to count matrix rows
anno <- anno[match(rownames(counts), anno$gene_id), ]

## 3. Filter to well annotated protein coding genes

# Predicted genes: Gm* and anything ending with "Rik" optionally followed by digits
is_predicted <- grepl("^Gm[0-9]+", anno$gene_name) |
                grepl("Rik$", anno$gene_name)
is_predicted[is.na(is_predicted)] <- FALSE

keep <- (anno$biotype == "protein_coding") & !is_predicted
keep[is.na(keep)] <- FALSE

counts_pc <- counts[keep, ]
symbols   <- anno$gene_name[keep]

## 4. Aggregate to unique gene symbols

# Remove rows with missing symbols
valid <- !is.na(symbols) & symbols != ""
counts_pc <- counts_pc[valid, ]
symbols   <- symbols[valid]

# Aggregate counts for genes that share the same symbol
counts_symbol <- rowsum(counts_pc, group = symbols)

## 5. Remove unexpressed genes (all zero counts across samples)

expr_keep <- rowSums(counts_symbol) > 0
counts_symbol <- counts_symbol[expr_keep, ]

## 6. Round counts and save final matrix

counts_symbol_int <- round(counts_symbol)

write.table(
  counts_symbol_int,
  file = "gene_counts.tsv",
  sep = "\t",
  quote = FALSE,
  col.names = NA
)
