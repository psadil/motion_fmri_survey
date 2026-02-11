get_ukb <- function(src, exclusion) {
  arrow::open_dataset(src) |>
    dplyr::filter(t > 0) |>
    dplyr::mutate(
      time = t * 0.735,
      ped = "AP",
      scan = 1L,
    ) |>
    dplyr::collect() |>
    add_run() |>
    do_casting() |>
    dplyr::anti_join(exclusion, by = dplyr::join_by(sub, ses, task))
}

get_ukb_exclusion <- function(src) {
  readr::read_tsv(src, col_types = "cccc") |>
    dplyr::rename(sub = subject_id)
}


get_ukb_design <- function() {
  # readr::read_tsv(
  #   "data/ukb/design.mat",
  #   skip = 5,
  #   col_names = FALSE,
  #   col_select = c("X1", "X2", "X3", "X4"),
  #   col_types = "dddd") |>
  #   dplyr::mutate(
  #     dplyr::across(
  #       tidyselect::everything(),
  #       \(x) scale(x) * .1),
  #     t = 1:dplyr::n()
  #     )

  # https://git.fmrib.ox.ac.uk/falmagro/UK_biobank_pipeline_v_1/-/tree/master/bb_data
  # tibble::tibble(
  #   faces = c(21, 50, 121, 150, 221) / .735,
  #   shapes = c(0, 71, 100, 171, 200) / .735,
  #   duration = 21 / .735
  # ) |>
  #   tidyr::pivot_longer(
  #     c(faces, shapes),
  #     values_to = "onset",
  #     names_to = "Block")
  #
  readr::read_tsv(
    "https://git.fmrib.ox.ac.uk/falmagro/UK_biobank_pipeline_v_1/-/raw/master/bb_data/task-hariri_events.tsv?inline=false",
    show_col_types = FALSE
  ) |>
    dplyr::select(onset, duration, type = Stimulus)
}


get_ukb_responses <- function(src) {
  txt <- rprime::read_eprime(src) |>
    rprime::FrameList()

  txt |>
    rprime::to_data_frame() |>
    tibble::as_tibble() |>
    dplyr::filter(Running == "SessionSelectionList") |>
    dplyr::select(
      Procedure,
      Sample,
      RunTrialNumber,
      starts_with("Experimenter")
    ) |>
    dplyr::mutate(
      onset = as.numeric(ExperimenterWindow.OnsetTime),
      onset = onset - min(onset, na.rm = TRUE),
      onset = onset / 1000,
      task = "faces/shapes"
    ) |>
    dplyr::filter(Procedure == "TrialsPROC") |>
    dplyr::select(onset, task)
}

get_ukb_demographics <- function(bulk, ukb_source) {
  subses <- duckplyr::read_parquet_duckdb(ukb_source) |>
    dplyr::distinct(sub, ses) |>
    dplyr::mutate(dplyr::across(
      tidyselect::all_of(c("sub", "ses")),
      as.integer
    )) |>
    dplyr::collect()

  sites <- duckplyr::read_parquet_duckdb(bulk) |>
    dplyr::select(
      sub = eid,
      site_2 = f.54.2.0,
      site_3 = f.54.3.0
    ) |>
    tidyr::pivot_longer(
      c(tidyselect::ends_with("_2"), tidyselect::ends_with("_3")),
      names_to = "ses",
      values_to = "site",
      names_prefix = "site_"
    ) |>
    dplyr::mutate(ses = as.integer(ses))

  ages <- duckplyr::read_parquet_duckdb(bulk) |>
    dplyr::select(
      sub = eid,
      age_2 = f.21003.2.0,
      age_3 = f.21003.3.0
    ) |>
    tidyr::pivot_longer(
      c(tidyselect::ends_with("_2"), tidyselect::ends_with("_3")),
      names_to = "ses",
      values_to = "age",
      names_prefix = "age_",
      values_ptypes = numeric()
    ) |>
    dplyr::mutate(ses = as.integer(ses))

  duckplyr::read_parquet_duckdb(bulk) |>
    dplyr::select(
      sub = eid,
      sex = f.31.0.0,
      tidyselect::starts_with("f.21000"),
      bmi_2 = f.21001.2.0,
      bmi_3 = f.21001.3.0
    ) |>
    dplyr::rename(ethnicity = f.21000.0.0) |>
    dplyr::mutate(
      ethnicity = dplyr::if_else(is.na(ethnicity), f.21000.1.0, ethnicity)
    ) |>
    dplyr::mutate(
      ethnicity = dplyr::if_else(is.na(ethnicity), f.21000.2.0, ethnicity)
    ) |>
    dplyr::select(-tidyselect::starts_with("f.21000")) |>
    tidyr::pivot_longer(
      c(tidyselect::ends_with("_2"), tidyselect::ends_with("_3")),
      names_to = "ses",
      values_to = "bmi",
      names_prefix = "bmi_",
      values_ptypes = numeric()
    ) |>
    dplyr::mutate(ses = as.integer(ses)) |>
    dplyr::semi_join(subses, by = dplyr::join_by(sub, ses)) |>
    dplyr::left_join(sites, by = dplyr::join_by(sub, ses)) |>
    dplyr::left_join(ages, by = dplyr::join_by(sub, ses)) |>
    do_casting()
}
