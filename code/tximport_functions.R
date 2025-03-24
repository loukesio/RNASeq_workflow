library(tximport)
import_kallisto_to_gene_level <- function(kallisto_dir, tx2gene, add_prefix = TRUE, file_pattern = "abundance.h5", save_path = NULL) {
  if (is(tx2gene, "TxDb") || (is.character(tx2gene) && file.exists(tx2gene))) tx2gene <- create_tx2gene(tx2gene, add_prefix)
  dirs <- list.files(kallisto_dir, full.names = TRUE)
  files <- file.path(dirs, file_pattern)
  names(files) <- basename(dirs)
  files <- files[file.exists(files)]
  txi <- tximport(files, type = "kallisto", tx2gene = tx2gene, countsFromAbundance = "lengthScaledTPM")
  if (!is.null(save_path)) {
    dir.create(dirname(save_path), recursive = TRUE, showWarnings = FALSE)
    save(txi, file = save_path)
  }
  return(txi)
}