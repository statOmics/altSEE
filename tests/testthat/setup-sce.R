# Minimal SCE for altSEE tests:
#   Main experiment: 4 proteins x 6 samples, with assay "proteins",
#   rowData columns Proteins (character), PValue, logFC, AveExpr (numeric),
#   reducedDim PCA.
#   altExp "peptides": 10 peptides x 6 samples, assay "peptides_norm",
#   rowData columns Proteins (character), PValue, logFC, AveExpr (numeric),
#   reducedDim PCA.
#   colData: Treatment (factor), SampleGroup (character).

library(SingleCellExperiment)
library(SummarizedExperiment)

set.seed(42)
n_prot  <- 4L
n_pep   <- 10L
n_samp  <- 6L

prot_names <- paste0("PROT", seq_len(n_prot))
pep_names  <- paste0("PEP",  seq_len(n_pep))
samp_names <- paste0("S",    seq_len(n_samp))

# peptide -> protein mapping: first 3 peptides map to PROT1, next 3 to PROT2,
# next 2 to PROT3, last 2 to PROT4
pep_protein <- c(rep("PROT1", 3), rep("PROT2", 3), rep("PROT3", 2), rep("PROT4", 2))

prot_mat <- matrix(rnorm(n_prot * n_samp), nrow = n_prot, ncol = n_samp,
                   dimnames = list(prot_names, samp_names))
pep_mat  <- matrix(rnorm(n_pep  * n_samp), nrow = n_pep,  ncol = n_samp,
                   dimnames = list(pep_names,  samp_names))

prot_rd <- DataFrame(
    Proteins  = prot_names,
    PValue    = runif(n_prot, 0.001, 0.5),
    logFC     = rnorm(n_prot),
    AveExpr   = rnorm(n_prot),
    row.names = prot_names
)

pep_rd <- DataFrame(
    Proteins  = pep_protein,
    PValue    = runif(n_pep, 0.001, 0.5),
    logFC     = rnorm(n_pep),
    AveExpr   = rnorm(n_pep),
    row.names = pep_names
)

col_df <- DataFrame(
    Treatment   = factor(rep(c("A", "B"), each = n_samp / 2)),
    SampleGroup = rep(c("ctrl", "trt"),  each = n_samp / 2),
    row.names   = samp_names
)

pep_sce <- SingleCellExperiment(
    assays   = list(peptides_norm = pep_mat),
    rowData  = pep_rd,
    colData  = col_df
)
reducedDim(pep_sce, "PCA") <- matrix(rnorm(n_samp * 4), nrow = n_samp,
    dimnames = list(samp_names, paste0("PC", 1:4)))

sce <- SingleCellExperiment(
    assays   = list(proteins = prot_mat),
    rowData  = prot_rd,
    colData  = col_df
)
reducedDim(sce, "PCA") <- matrix(rnorm(n_samp * 4), nrow = n_samp,
    dimnames = list(samp_names, paste0("PC", 1:4)))
altExp(sce, "peptides") <- pep_sce
