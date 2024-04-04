get_ukb <- function(sources){
  
  # batch reading. Otherwise errors about too many open files
  batches <- tibble::tibble(
    src = sources
  ) |>
    dplyr::mutate(
      batch = 1:dplyr::n() %% 10
    ) |>
    dplyr::group_nest(batch)
  
  purrr::map(
    batches$data,
    ~readr::read_delim(
      .x$src,
      col_names = c("rot_x", "rot_y", "rot_z", "trans_x", "trans_y", "trans_z", "t"),
      col_select = c("rot_x", "rot_y", "rot_z", "trans_x", "trans_y", "trans_z", "t"),
      col_types = "ddddddi",
      delim = "  ",
      id = "path",
      trim_ws = TRUE) 
  ) |>
    dplyr::bind_rows() |>
    dplyr::mutate(
      sub = stringr::str_extract(path, "[[:digit:]]{7}") |>
        as.integer(),
      ses = stringr::str_extract(path, "(?<=ses-)[[:digit:]]") |>
        as.integer(),
      .by = c(path)
    ) |>
    dplyr::select(-path)
}

get_ukb_reg <- function(ukb, include_ukb){
  ukb |>
    dplyr::mutate(ntr = dplyr::n(), .by = c(sub, ses)) |>
    dplyr::filter(ntr == 332) |>
    dplyr::select(-ntr) |>
    dplyr::semi_join(include_ukb, by = dplyr::join_by(sub))
}

get_ukb_design <- function(mat){
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
  readr::read_tsv("https://git.fmrib.ox.ac.uk/falmagro/UK_biobank_pipeline_v_1/-/raw/master/bb_data/task-hariri_events.tsv?inline=false") |>
    dplyr::select(onset, duration, Block=Stimulus) |>
    mutate(across(c(onset, duration), \(x) x / 0.735))
}

get_include_ukb <- function(path){
  polars::pl$scan_parquet(
    path
  )$select(
    "eid"
  )$rename(
    sub="eid"
  )$collect()$to_data_frame() |>
    tibble::as_tibble()
}

get_ukb_responses <- function(src="data/1000513_25748_2_0.txt"){
  txt <- rprime::read_eprime(src) |>
    rprime::FrameList()
  
  txt |>
    rprime::to_data_frame() |> 
    tibble::as_tibble() |> 
    dplyr::filter(Running=="SessionSelectionList") |> 
    dplyr::select(Procedure, Sample, RunTrialNumber, starts_with("Experimenter")) |> 
    dplyr::mutate(
      onset = as.numeric(ExperimenterWindow.OnsetTime), 
      onset = onset - min(onset, na.rm = TRUE),
      onset = onset / 1000 / 0.735) |> 
    dplyr::filter(Procedure=="TrialsPROC") |>
    dplyr::select(onset)
}
