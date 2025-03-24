validate_exp_design <- function(exp_design, required_cols = c("sample_name", "treatment", "species")) {
  if (!all(required_cols %in% colnames(exp_design))) stop("Missing required columns")
  if(any(duplicated(exp_design$sample_name))) stop("Duplicate sample names detected")
  message("Design validation successful")
  invisible(exp_design)
}