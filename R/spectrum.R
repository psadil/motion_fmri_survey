rescale <- function(.data, ...) {
  .data |>
    dplyr::mutate(pxx = 10 * log10(pxx)) |>
    na.omit() |>
    dplyr::group_nest(...) |>
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
    dplyr::select(-lower, -upper)
}

.avg_by_run <- function(by_run) {
  by_run |>
    dplyr::filter(!filtered) |>
    dplyr::summarise(avg = median(loc), .by = c(sub, dataset))
}

.clean <- function(.data, N = 100) {
  .data |> dplyr::semi_join(dplyr::count(.data, freq) |> dplyr::filter(n > N))
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

get_abcd_spectrum <- function(src, by_run) {
  to_keep <- by_run |>
    dplyr::filter(dataset == "abcd") |>
    dplyr::distinct(sub, task, ses, scan) |>
    dplyr::mutate(scan = as.integer(scan))

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
      dplyr::rename(scan = run) |>
      do_casting() |>
      dplyr::collect() |>
      dplyr::inner_join(to_keep, by = dplyr::join_by(sub, ses, task, scan))
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

bind_spectrum <- function(spectrums) {
  dplyr::bind_rows(spectrums, .id = "dataset") |> dplyr::select(-ped, -run)
}

rescale_spectrum <- function(spectrums, by_run) {
  spectrums |>
    rescale(param, dataset, sub, ses, task, scan) |>
    dplyr::inner_join(.avg_by_run(by_run), by = dplyr::join_by(sub, dataset))
}

get_ukb_spectrum <- function(src, ukb) {
  duckplyr::read_parquet_duckdb(src, prudence = "lavish") |>
    dplyr::mutate(ped = "AP", scan = 1, run = 1) |>
    do_casting() |>
    dplyr::inner_join(
      dplyr::distinct(ukb, sub, task, run, ped, ses, scan),
      by = dplyr::join_by(sub, ses, task, ped, run, scan)
    )
}

get_limits <- function(
  d,
  lower_exclude = 0.15,
  upper_exclude = 0.6,
  trans_dir = "trans_y"
) {
  d |>
    dplyr::filter(
      dplyr::between(freq, lower_exclude, upper_exclude),
      param == trans_dir
    ) |>
    dplyr::slice_max(order_by = pxx, n = 1, by = c(sub, ses, task, scan)) |>
    dplyr::summarise(
      lower = median(freq),
      upper = quantile(freq, 0.75),
      .by = c(task)
    ) |>
    dplyr::summarise(lower = mean(lower), upper = mean(upper))
}

prep_tbl_freq_limits <- function(
  hcpa_spectrum,
  hcpd_spectrum,
  hcpya_spectrum,
  spacetop_spectrum,
  ukb_spectrum,
  abcd_spectrum,
  demographics
) {
  hcpd_young <- hcpd_spectrum |>
    dplyr::left_join(
      demographics |>
        dplyr::filter(dataset == "hcpd") |>
        dplyr::distinct(sub, ses, age)
    ) |>
    dplyr::filter(age < 8) |>
    get_limits() |>
    dplyr::mutate(dataset = "hcpd", ses = "6-7")

  hcpd_old <- hcpd_spectrum |>
    dplyr::left_join(
      demographics |>
        dplyr::filter(dataset == "hcpd") |>
        dplyr::distinct(sub, ses, age)
    ) |>
    dplyr::filter(age >= 8) |>
    get_limits() |>
    dplyr::mutate(dataset = "hcpd", ses = "8+")

  abcd <- purrr::map(
    list(
      Year2 = dplyr::filter(abcd_spectrum, ses == "Year2"),
      Year4 = dplyr::filter(abcd_spectrum, ses == "Year4"),
      Baseline = dplyr::filter(abcd_spectrum, ses == "Baseline")
    ),
    get_limits
  ) |>
    dplyr::bind_rows(.id = "ses") |>
    dplyr::mutate(dataset = "abcd")

  purrr::map2(
    list(
      hcpya = hcpya_spectrum,
      hcpa = hcpa_spectrum,
      spacetop = spacetop_spectrum,
      ukb = ukb_spectrum
    ),
    c("trans_x", "trans_y", "trans_y", "trans_y"),
    ~ get_limits(.x, trans_dir = .y)
  ) |>
    dplyr::bind_rows(.id = "dataset") |>
    dplyr::bind_rows(hcpd_young, hcpd_old, abcd)
}
