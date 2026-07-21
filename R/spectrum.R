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


get_hcpa_spectrum <- function(src, hcpa, by_run) {
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
    ) |>
    dplyr::select(-ped, -run) |>
    rescale(param, sub, ses, task, scan) |>
    dplyr::mutate(dataset = "hcpa") |>
    dplyr::inner_join(.avg_by_run(by_run), by = dplyr::join_by(sub, dataset))
}

get_hcpd_spectrum <- function(src, hcpd, by_run) {
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
    ) |>
    dplyr::select(-ped, -run) |>
    dplyr::filter(scan == 1) |>
    rescale(param, sub, ses, task, scan) |>
    dplyr::mutate(dataset = "hcpd") |>
    dplyr::inner_join(.avg_by_run(by_run), by = dplyr::join_by(sub, dataset))
}

get_hcpya_spectrum <- function(src, hcpya, by_run) {
  rest_exclude <- duckplyr::read_parquet_duckdb(src, prudence = "lavish") |>
    dplyr::filter(ped == "RL") |>
    dplyr::count(task, freq) |>
    dplyr::filter(n < 10) |>
    dplyr::mutate(task = "rest")

  duckplyr::read_parquet_duckdb(src, prudence = "lavish") |>
    dplyr::filter(ped == "RL") |>
    dplyr::mutate(
      task = dplyr::case_when(
        stringr::str_detect(task, "rest") ~ "rest",
        .default = task
      )
    ) |>
    do_casting() |>
    dplyr::semi_join(
      dplyr::distinct(hcpya, sub, task, ped),
      by = dplyr::join_by(sub, task, ped)
    ) |>
    dplyr::anti_join(rest_exclude, by = dplyr::join_by(task, freq)) |>
    dplyr::select(-ped) |>
    rescale(param, sub, ses, task) |>
    dplyr::mutate(dataset = "hcpya") |>
    dplyr::inner_join(.avg_by_run(by_run), by = dplyr::join_by(sub, dataset)) |>
    dplyr::mutate(scan = 1)
}

get_abcd_spectrum <- function(src, by_run) {
  # for simplicity, keep only participants with the most typical
  # frequencies (different scan lengths result in diff freqs)

  to_keep <- by_run |>
    dplyr::filter(dataset == "abcd") |>
    dplyr::distinct(sub, task, ses, scan) |>
    dplyr::mutate(scan = as.integer(scan))

  freqs <- duckplyr::read_parquet_duckdb(src, prudence = "lavish") |>
    dplyr::filter(run == 1) |>
    dplyr::count(ses, task, freq) |>
    dplyr::filter(n > 10)

  duckplyr::read_parquet_duckdb(src, prudence = "lavish") |>
    dplyr::rename(scan = run) |>
    dplyr::filter(scan == 1) |>
    dplyr::semi_join(to_keep, by = dplyr::join_by(sub, ses, task, scan)) |>
    dplyr::semi_join(freqs, by = dplyr::join_by(ses, task, freq)) |>
    dplyr::collect() |>
    rescale(param, sub, ses, task, scan) |>
    dplyr::mutate(dataset = "abcd") |>
    dplyr::inner_join(.avg_by_run(by_run), by = dplyr::join_by(sub, dataset))
}

get_spacetop_spectrum <- function(src, spacetop, by_run) {
  freqs <- duckplyr::read_parquet_duckdb(src, prudence = "lavish") |>
    dplyr::mutate(scan = as.integer(run)) |>
    dplyr::filter(scan == 1) |>
    dplyr::count(task, freq)

  narratives <- freqs |>
    dplyr::filter(task == "narratives") |>
    dplyr::filter(n == 624)

  fractional <- freqs |>
    dplyr::filter(task == "fractional") |>
    dplyr::filter(n == 588)

  alignvideo <- freqs |>
    dplyr::filter(task == "alignvideo") |>
    dplyr::filter(n == 618)

  # shortvideo has two comparable groups of participants
  # shortvideo <- freqs |>
  #   dplyr::semi_join(freq_short, by = dplyr::join_by(task, freq))
  shortvideo <- duckplyr::read_parquet_duckdb(src, prudence = "lavish") |>
    dplyr::mutate(scan = as.integer(run)) |>
    dplyr::filter(scan == 1, task == "shortvideo") |>
    dplyr::filter(as.integer(sub) < 31) |>
    dplyr::distinct(task, sub)

  duckplyr::read_parquet_duckdb(src, prudence = "lavish") |>
    dplyr::filter(freq > 0) |>
    dplyr::mutate(ses = as.integer(ses), scan = run) |> # need to strip leading 0 before do_casting
    do_casting() |>
    dplyr::filter(scan == 1) |>
    dplyr::anti_join(shortvideo, by = dplyr::join_by(task, sub)) |>
    dplyr::semi_join(
      dplyr::bind_rows(
        dplyr::filter(
          freqs,
          !(task %in% c("alignvideo", "fractional", "narratives"))
        ),
        alignvideo,
        fractional,
        narratives
      ),
      by = dplyr::join_by(task, freq)
    ) |>
    dplyr::inner_join(
      dplyr::distinct(spacetop, sub, task, run, ses, scan),
      by = dplyr::join_by(sub, ses, task, run, scan)
    ) |>
    dplyr::select(-run) |>
    rescale(param, sub, ses, task, scan) |>
    dplyr::mutate(dataset = "spacetop") |>
    dplyr::inner_join(.avg_by_run(by_run), by = dplyr::join_by(sub, dataset))
}

bind_spectrum <- function(spectrums, by_run) {
  dplyr::bind_rows(spectrums, .id = "dataset") |>
    dplyr::inner_join(.avg_by_run(by_run), by = dplyr::join_by(sub, dataset))
}

rescale_spectrum <- function(spectrums, by_run) {
  spectrums |>
    rescale(param, dataset, sub, ses, task, scan) |>
    dplyr::inner_join(.avg_by_run(by_run), by = dplyr::join_by(sub, dataset))
}

get_ukb_spectrum <- function(src, ukb, by_run) {
  duckplyr::read_parquet_duckdb(src, prudence = "lavish") |>
    dplyr::mutate(scan = 1) |>
    do_casting() |>
    dplyr::inner_join(
      dplyr::distinct(ukb, sub, task, ses, scan),
      by = dplyr::join_by(sub, ses, task, scan)
    ) |>
    dplyr::filter(ses == "3") |>
    rescale(param, sub, ses, task, scan) |>
    dplyr::mutate(dataset = "ukb") |>
    dplyr::inner_join(.avg_by_run(by_run), by = dplyr::join_by(sub, dataset))
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
      Baseline = dplyr::filter(abcd_spectrum, ses == "Baseline"),
      Year6 = dplyr::filter(abcd_spectrum, ses == "Year6")
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
