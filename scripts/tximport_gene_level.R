# ------------------------------------------------------------
# Gene-level summarisation of Kallisto output using tximport
# ------------------------------------------------------------

library(tximport)

# Project folders
kallisto_dir <- "data/quant"
results_dir <- "results"

# ------------------------------------------------------------
# 1. Transcript-to-gene mapping
# ------------------------------------------------------------

t2g <- read.csv(
    "t2g.txt",
    sep = "\t",
    header = FALSE
)

# tximport only needs transcript ID and gene ID
tx2gene <- t2g[, c(1, 2)]

colnames(tx2gene) <- c(
    "transcript_id",
    "gene_id"
)

# Remove Ensembl gene version numbers if present
tx2gene$gene_id <- sub(
    "\\.[0-9]+$",
    "",
    tx2gene$gene_id
)

# ------------------------------------------------------------
# 2. Find all Kallisto samples
# ------------------------------------------------------------

sample_dirs <- list.dirs(
    kallisto_dir,
    recursive = FALSE,
    full.names = TRUE
)

sample_names <- basename(sample_dirs)

# Path to abundance.h5 for each sample
paths <- file.path(
    sample_dirs,
    "abundance.h5"
)

names(paths) <- sample_names

print(paths)

# Check that every file exists
stopifnot(all(file.exists(paths)))

# ------------------------------------------------------------
# 3. Run tximport
# ------------------------------------------------------------

txi <- tximport(
    paths,
    type = "kallisto",
    tx2gene = tx2gene,
    countsFromAbundance = "lengthScaledTPM"
)

# ------------------------------------------------------------
# 4. Extract gene-level matrices
# ------------------------------------------------------------

gene_tpm_all <- txi$abundance
gene_counts_all <- txi$counts

print(dim(gene_tpm_all))
print(dim(gene_counts_all))

# ------------------------------------------------------------
# 5. Save results
# ------------------------------------------------------------

write.csv(
    gene_tpm_all,
    file.path(
        results_dir,
        "gene_tpm_all_tximport.csv"
    ),
    quote = FALSE
)

write.csv(
    gene_counts_all,
    file.path(
        results_dir,
        "gene_counts_all_tximport.csv"
    ),
    quote = FALSE
)

cat("tximport finished successfully\n")