# R/standardize_input.R

#' Standardize column names using the provided meta information file
#' @param data Input dataset containing merged metabolic and transcriptomic data
#' @param info Metadata for column standardization
#' @export
#'
standardize_input <- function(data, info){
  m = match(colnames(data), info$Variable)
  m = m[!is.na(m)]
  colnames(data)[colnames(data) %in% info$Variable] = info$Colname[m]
  return(data)
}
