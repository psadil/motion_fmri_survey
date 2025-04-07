get_mriqc_bold <- function() {
  arrow::open_dataset("data/bold") |>
    dplyr::collect()
}
