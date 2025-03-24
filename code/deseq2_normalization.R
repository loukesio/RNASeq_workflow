library(DESeq2)
run_deseq2_vst <- function(txi, sample_table, design_formula = ~condition, save_path = NULL) {
  dds <- DESeqDataSetFromMatrix(round(txi$counts), colData = sample_table, design = design_formula)
  dds <- estimateSizeFactors(dds)
  vst_counts <- vst(dds, blind = FALSE)
  ml_data <- assay(vst_counts)
  if (!is.null(save_path)) {
    dir.create(dirname(save_path), recursive = TRUE, showWarnings = FALSE)
    save(ml_data, file = save_path)
  }
  return(ml_data)
}