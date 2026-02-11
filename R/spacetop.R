get_spacetop_demographics <- function(src) {
  readr::read_tsv(src, col_select = c(-group), show_col_types = FALSE) |>
    dplyr::rename(sub = participant_id) |>
    dplyr::mutate(
      sex = dplyr::replace_values(
        sex,
        "M" ~ "Male",
        "F" ~ "Female"
      ),
      sub = stringr::str_extract(sub, "(?<=sub-)[[:digit:]]+")
    ) |>
    do_casting()
}

get_spacetop <- function(file, exclusion) {
  duckplyr::read_parquet_duckdb(file, prudence = "lavish") |>
    dplyr::mutate(ses = as.integer(ses)) |> # need to strip leading 0 before do_casting
    do_casting() |>
    dplyr::anti_join(exclusion, by = dplyr::join_by(sub, ses, task, run)) |>
    dplyr::mutate(time = t * 0.46, ped = "AP", scan = run) |>
    dplyr::collect()
}

get_spacetop_exclusion <- function(
  src = "data/exclusion/Spacetop_exclusions.rds"
) {
  exclusion <- readRDS(src)
  exclusion$Spacetop_exclusion_table |>
    tibble::as_tibble() |>
    dplyr::select(
      sub = Excluded_sub,
      ses = Excluded_ses,
      task = Excluded_task,
      run = Excluded_run
    ) |>
    do_casting()
}
