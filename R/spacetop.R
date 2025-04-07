get_spacetop_demographics <- function() {
  readr::read_tsv(here::here("data/participants.tsv"), col_select = c(-group)) |>
    dplyr::rename(sub = participant_id) |>
    dplyr::mutate(
      sex = dplyr::case_match(
        sex,
        "M" ~ "Male",
        "F" ~ "Female"
      )
    ) |>
    dplyr::mutate(sub = as.character(sub))
}

get_spacetop <- function(file, exclusion) {
  arrow::open_dataset(file) |>
    dplyr::collect() |>
    dplyr::anti_join(exclusion, by = dplyr::join_by(sub, ses, task, run)) |>
    dplyr::mutate(
      time = t * 0.46,
      ped = "AP",
      ses = as.integer(ses),
      scan = as.integer(run),
      run = as.integer(run)
    ) |>
    do_casting()
}

get_spacetop_exclusion <- function(src = "data/exclusion/Spacetop_exclusions.rds") {
  exclusion <- readRDS(src)
  exclusion$Spacetop_exclusion_table |>
    tibble::as_tibble() |>
    dplyr::select(sub = Excluded_sub, ses = Excluded_ses, task = Excluded_task, run = Excluded_run) |>
    dplyr::mutate(run = as.character(run))
}
