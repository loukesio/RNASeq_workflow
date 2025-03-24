library(GenomicFeatures)
create_txdb <- function(gff_file, organism, save_path = NULL) {
  txdb <- makeTxDbFromGFF(file = gff_file, dataSource = paste("GFF file for", organism), organism = organism)
  if (!is.null(save_path)) {
    dir.create(dirname(save_path), recursive = TRUE, showWarnings = FALSE)
    saveDb(txdb, save_path)
  }
  return(txdb)
}

load_or_create_txdb <- function(gff_file, organism, db_path, force_recreate = FALSE) {
  if (file.exists(db_path) && !force_recreate) txdb <- loadDb(db_path)
  else txdb <- create_txdb(gff_file, organism, db_path)
  return(txdb)
}

create_tx2gene <- function(txdb, add_prefix = TRUE) {
  if (is.character(txdb) && file.exists(txdb)) txdb <- loadDb(txdb)
  keys_db <- keys(txdb, keytype = "GENEID")
  tx2gene <- select(txdb, keys = keys_db, keytype = "GENEID", columns = "TXNAME")
  colnames(tx2gene) <- c("GENEID", "target_id")
  if (add_prefix) tx2gene$target_id <- paste0("transcript:", tx2gene$target_id)
  return(tx2gene)
}