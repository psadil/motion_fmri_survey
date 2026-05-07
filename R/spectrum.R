rescale <- function(.data, ...) {
  .data |>
    dplyr::mutate(
      avg = factor(avg, ordered = TRUE),
      pxx = 10 * log10(pxx)
    ) |>
    na.omit() |>
    dplyr::group_nest(..., avg) |>
    dplyr::mutate(
      data = purrr::map(
        data,
        ~ tibble::tibble(
          pxx = lm(scale(pxx) ~ pxx, data = .x) |> predict(),
          freq = .x$freq,
        )
      ),
      lower = purrr::map_dbl(data, ~ quantile(.x$pxx, 0.01)),
      upper = purrr::map_dbl(data, ~ quantile(.x$pxx, 0.96))
    ) |>
    tidyr::unnest(data) |>
    dplyr::mutate(
      pxx = (pxx - min(pxx)) / (max(pxx) - min(pxx)),
      pxx = dplyr::if_else(pxx < lower, lower, pxx),
      pxx = dplyr::if_else(pxx > upper, upper, pxx),
      .by = c(...)
    ) |>
    dplyr::mutate(avg = as.numeric(as.character(avg))) |>
    dplyr::select(-lower, -upper)
}

.avg_by_run <- function(by_run) {
  by_run |>
    dplyr::filter(!filtered) |>
    dplyr::summarise(
      avg = median(loc),
      .by = c(sub, dataset)
    )
}

.clean <- function(.data, N = 100) {
  .data |>
    dplyr::semi_join(dplyr::count(.data, freq) |> dplyr::filter(n > N))
}


get_hcpa_spectrum <- function(src, hcpa) {
  duckplyr::read_parquet_duckdb(src, prudence = "lavish") |>
    add_run() |>
    dplyr::mutate(
      task = stringr::str_to_lower(task),
      run = dplyr::recode_values(
        task,
        "rest1" ~ 1L,
        "rest2" ~ 2L,
        default = run
      )
    ) |>
    do_casting() |>
    dplyr::inner_join(
      dplyr::distinct(hcpa, sub, task, run, ped, ses, scan),
      by = dplyr::join_by(sub, ses, task, ped, run)
    )
}

get_hcpd_spectrum <- function(src, hcpd) {
  duckplyr::read_parquet_duckdb(src, prudence = "lavish") |>
    add_run() |>
    dplyr::mutate(
      task = stringr::str_to_lower(task),
      run = dplyr::recode_values(
        task,
        "rest1" ~ 1L,
        "rest2" ~ 2L,
        default = run
      ),
      task = dplyr::case_when(
        task == "rest1a" ~ "resta",
        task == "rest1b" ~ "restb",
        task == "rest2a" ~ "resta",
        task == "rest2b" ~ "restb",
        stringr::str_detect(task, "rest") ~ "rest",
        .default = task
      )
    ) |>
    do_casting() |>
    dplyr::inner_join(
      dplyr::distinct(hcpd, sub, task, run, ped, ses, scan),
      by = dplyr::join_by(sub, ses, task, ped, run)
    )
}

get_hcpya_spectrum <- function(src, hcpya) {
  duckplyr::read_parquet_duckdb(src, prudence = "lavish") |>
    add_run() |>
    dplyr::mutate(
      task = stringr::str_to_lower(task),
      run = dplyr::recode_values(
        task,
        "rest1" ~ 1L,
        "rest2" ~ 2L,
        default = run
      ),
      task = dplyr::case_when(
        task == "rest1a" ~ "resta",
        task == "rest1b" ~ "restb",
        task == "rest2a" ~ "resta",
        task == "rest2b" ~ "restb",
        stringr::str_detect(task, "rest") ~ "rest",
        .default = task
      )
    ) |>
    do_casting() |>
    dplyr::inner_join(
      dplyr::distinct(hcpya, sub, task, run, ped, ses, scan),
      by = dplyr::join_by(sub, ses, task, ped, run)
    )
}

get_abcd_spectrum <- function(src, abcd) {
  to_keep <- dplyr::distinct(abcd, sub, task, run, ped, ses, scan)

  unique_ses_src <- duckplyr::read_parquet_duckdb(src) |>
    dplyr::distinct(ses, run) |>
    dplyr::collect()
  out <- vector("list", nrow(unique_ses_src))
  for (i in seq_len(nrow(unique_ses_src))) {
    out[[i]] <- duckplyr::read_parquet_duckdb(src, prudence = "lavish") |>
      dplyr::filter(
        ses == unique_ses_src$ses[[i]],
        run == unique_ses_src$run[[i]]
      ) |>
      convert_abcd_ses() |>
      do_casting() |>
      dplyr::inner_join(to_keep, by = dplyr::join_by(sub, ses, task, run))
  }
  dplyr::bind_rows(out)
}

get_spacetop_spectrum <- function(src, spacetop) {
  duckplyr::read_parquet_duckdb(src, prudence = "lavish") |>
    dplyr::mutate(ses = as.integer(ses), ped = "AP", scan = run) |> # need to strip leading 0 before do_casting
    do_casting() |>
    dplyr::inner_join(
      dplyr::distinct(spacetop, sub, task, run, ped, ses, scan),
      by = dplyr::join_by(sub, ses, task, ped, run, scan)
    )
}

bind_spectrum <- function(datasets, by_run) {
  dplyr::bind_rows(datasets, .id = "dataset") |>
    dplyr::inner_join(
      .avg_by_run(by_run),
      by = dplyr::join_by(sub, dataset)
    ) |>
    rescale(task, ses, run, sub, param, scan, dataset)
}
