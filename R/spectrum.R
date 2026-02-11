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


get_hcpa_spectrum <- function(src, by_run) {
  duckplyr::read_parquet_duckdb(src, prudence = "lavish") |>
    do_casting() |>
    dplyr::left_join(
      .avg_by_run(dplyr::filter(by_run, dataset == "hcpa")),
      by = dplyr::join_by(sub)
    ) |>
    dplyr::filter(!(task == "REST1" & ped == "PA")) |>
    rescale(task, ped, sub, param)
}

get_hcpd_spectrum <- function(src, by_run, excluded) {
  duckplyr::read_parquet_duckdb(src, prudence = "lavish") |>
    do_casting() |>
    dplyr::left_join(
      .avg_by_run(dplyr::filter(by_run, dataset == "hcpd")),
      by = dplyr::join_by(sub)
    ) |>
    dplyr::filter(!(task %in% c("REST1a", "REST1b", "REST2a", "REST2b"))) |>
    dplyr::mutate(
      task = stringr::str_to_lower(task)
    ) |>
    dplyr::filter(!is.na(avg)) |>
    dplyr::anti_join(excluded) |>
    rescale(task, ped, sub, param)
}

get_hcpya_spectrum <- function(src, by_run, excluded) {
  duckplyr::read_parquet_duckdb(src, prudence = "lavish") |>
    do_casting() |>
    dplyr::left_join(
      .avg_by_run(dplyr::filter(by_run, dataset == "hcpya")),
      by = dplyr::join_by(sub)
    ) |>
    dplyr::filter(!is.na(avg)) |>
    dplyr::anti_join(excluded, by = dplyr::join_by(sub, task, ped)) |>
    rescale(task, ped, sub, param)
}

get_abcd_spectrum <- function(src, by_run, excluded) {
  by_run_abcd <- by_run |>
    dplyr::filter(dataset == "abcd", task == "rest", !filtered) |>
    dplyr::summarise(avg = median(loc), .by = c(sub, ses))

  unique_ses_src <- duckplyr::read_parquet_duckdb(src) |>
    dplyr::distinct(ses) |>
    purrr::pluck("ses")
  out <- vector("list", length(unique_ses_src))
  for (i in seq_along(unique_ses_src)) {
    out[[i]] <- duckplyr::read_parquet_duckdb(src, prudence = "lavish") |>
      dplyr::filter(ses == unique_ses_src[i]) |>
      convert_abcd_ses() |>
      do_casting() |>
      dplyr::anti_join(excluded, by = dplyr::join_by(sub, ses, task, run)) |>
      dplyr::left_join(by_run_abcd, by = dplyr::join_by(sub, ses)) |>
      dplyr::filter(!is.na(avg)) |>
      rescale(task, sub, ses, run, param)
  }
  dplyr::bind_rows(out)
}

get_spacetop_spectrum <- function(src, by_run, excluded) {
  duckplyr::read_parquet_duckdb(src, prudence = "lavish") |>
    do_casting() |>
    dplyr::left_join(
      .avg_by_run(dplyr::filter(by_run, dataset == "spacetop")),
      by = dplyr::join_by(sub)
    ) |>
    dplyr::filter(!is.na(avg)) |>
    dplyr::anti_join(excluded, by = dplyr::join_by(sub, ses, task, run)) |>
    rescale(task, ses, run, sub, param)
}
