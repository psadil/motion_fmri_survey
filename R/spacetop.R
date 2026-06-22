get_spacetop_demographics <- function(src) {
  readr::read_tsv(src, col_select = c(-group), show_col_types = FALSE) |>
    dplyr::rename(sub = participant_id) |>
    dplyr::mutate(sub = stringr::str_extract(sub, "(?<=sub-)[[:digit:]]+")) |>
    do_casting()
}

get_spacetop <- function(file) {
  duckplyr::read_parquet_duckdb(file, prudence = "lavish") |>
    dplyr::mutate(ses = as.integer(ses)) |> # need to strip leading 0 before do_casting
    do_casting() |>
    dplyr::mutate(time = t * 0.46, ped = "AP", scan = run) |>
    dplyr::collect() |>
    truncate_to_modal_lengths()
}
