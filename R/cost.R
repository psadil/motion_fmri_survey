get_params <- function(
  src = "https://thomasyeolab.github.io/OptimalScanTimeCalculator/budgetcalc.html"
) {
  html <- rvest::read_html_live(src) |>
    rvest::html_element("#phenotype-table") |>
    rvest::html_table() |>
    dplyr::rename_all(
      ~ stringr::str_remove_all(
        .x,
        "[[:space:]]*\\\n[[:space:]]*▼*\\.*[[:print:]]*"
      )
    ) |>
    dplyr::filter(stringr::str_detect(Dataset, "ABCD.*rest")) |>
    dplyr::summarise(K0 = mean(K0), K1 = mean(K1), K2 = mean(K2)) |>
    as.list()
  if (is.nan(html$K0)) {
    stop("Missing parameters")
  }
  html
}

accuracy <- function(n, t, k0, k1, k2) {
  k0 * sqrt(1 / (1 + k1 / n + k2 / (n * t)))
}

n_from_acc <- function(a, t = 30, k0, k1, k2) {
  (a^2 * (-k2) - a^2 * (k1) * t) / (t * (a^2 - k0^2))
}

cost <- function(
  n,
  t = 30,
  prop_lost = 0,
  cost_per_minute = 500 / 60,
  cost_per_participant = 500
) {
  # gemoetric series assumes we keep losing proportion
  (n / (1 - prop_lost)) * (cost_per_participant + cost_per_minute * t)
}

plot_cost <- function(by_run) {
  by_run |>
    dplyr::filter(
      task == "rest",
      ses == "baseline" | stringr::str_detect(dataset, "abcd", TRUE),
      ses == "2" | stringr::str_detect(dataset, "ukb", TRUE),
      !filtered,
      scan == "1"
    ) |>
    dplyr::mutate(N = dplyr::n(), .by = c(dataset, scan, ses, task)) |>
    dplyr::collect() |>
    tidyr::crossing(threshold = seq(0.1, 0.5, by = 0.01)) |>
    dplyr::mutate(kept = loc < threshold) |>
    dplyr::summarise(
      n = sum(kept),
      .by = c(dataset, task, scan, ses, threshold, N)
    ) |>
    dplyr::mutate(
      n_total = N / (n / N),
      t = dplyr::recode_values(
        dataset,
        c("abcd") ~ .8 * 382 / 60,
        c("hcpya") ~ .72 * 1200 / 60,
        "ukb" ~ .735 * 490 / 60,
        c("hcpa", "hcpd") ~ .72 * 478 / 60
      ),
      price = purrr::map2_dbl(n_total, t, cost) / cost(N, t),
      ideal = purrr::map_dbl(n_total, cost) / cost(N)
    ) |>
    ggplot2::ggplot(ggplot2::aes(x = threshold, y = price, color = dataset)) +
    ggplot2::geom_line() +
    ggplot2::scale_y_continuous(
      "Cost Factor Due\nto Exclusions",
      limits = c(1, NA),
      transform = "log2"
    ) +
    ggplot2::xlab("Threshold for\nAvg. Framewise Disp.") +
    ggplot2::scale_color_viridis_d(name = "Dataset", option = "turbo")
}

plot_cost_lost <- function(by_run) {
  by_run |>
    dplyr::filter(
      task == "rest",
      ses == "baseline" | stringr::str_detect(dataset, "abcd", TRUE),
      ses == "2" | stringr::str_detect(dataset, "ukb", TRUE),
      !filtered,
      scan == "1"
    ) |>
    dplyr::collect() |>
    dplyr::mutate(N = dplyr::n(), .by = c(dataset, scan, ses, task)) |>
    tidyr::crossing(threshold = seq(0.1, 0.5, by = 0.01)) |>
    dplyr::mutate(kept = loc < threshold) |>
    dplyr::summarise(
      n = sum(kept),
      .by = c(dataset, task, scan, ses, threshold, N)
    ) |>
    dplyr::mutate(
      t = dplyr::recode_values(
        dataset,
        c("abcd") ~ .8 * 382 / 60,
        c("hcpya") ~ .72 * 1200 / 60,
        "ukb" ~ .735 * 490 / 60,
        c("hcpa", "hcpd") ~ .72 * 478 / 60
      ),
      cost_prop = purrr::map2_dbl(n, t, cost) / purrr::map2_dbl(N, t, cost)
    ) |>
    ggplot2::ggplot(ggplot2::aes(
      x = threshold,
      y = cost_prop,
      color = dataset
    )) +
    ggplot2::geom_line() +
    ggplot2::scale_y_continuous(
      "Proportion of Budget\nProviding Usable Data",
      limits = c(0, 1)
    ) +
    ggplot2::xlab("Threshold for\nAvg. Framewise Disp.") +
    ggplot2::scale_color_viridis_d(name = "Dataset", option = "turbo")
}

plot_acc_decrease <- function(by_run, params) {
  by_run |>
    dplyr::filter(
      task == "rest",
      ses == "baseline" | stringr::str_detect(dataset, "abcd", TRUE),
      ses == "2" | stringr::str_detect(dataset, "ukb", TRUE),
      !filtered,
      scan == "1"
    ) |>
    dplyr::collect() |>
    dplyr::mutate(N = dplyr::n(), .by = c(dataset, scan, ses, task)) |>
    tidyr::crossing(threshold = seq(0.1, 0.5, by = 0.01)) |>
    dplyr::mutate(kept = loc < threshold) |>
    dplyr::summarise(
      n = sum(kept),
      .by = c(dataset, task, scan, ses, threshold, N)
    ) |>
    dplyr::mutate(
      n_total = N / (n / N),
      t = dplyr::recode_values(
        dataset,
        c("abcd") ~ .8 * 382 / 60,
        c("hcpya") ~ .72 * 1200 / 60,
        "ukb" ~ .735 * 490 / 60,
        c("hcpa", "hcpd") ~ .72 * 478 / 60
      ),
      acc = purrr::map2_dbl(
        n,
        t,
        accuracy,
        k0 = params$k0,
        k1 = params$k1,
        k2 = params$k2
      ) /
        accuracy(N, t, k0 = params$k0, k1 = params$k1, k2 = params$k2)
    ) |>
    ggplot2::ggplot(ggplot2::aes(x = threshold, y = acc, color = dataset)) +
    ggplot2::geom_line() +
    ggplot2::scale_y_continuous(
      "Model Performance\n(Scaled Correlation)",
      limits = c(0, 1)
    ) +
    ggplot2::xlab("Threshold for Avg\nFramewise Disp.") +
    ggplot2::scale_color_viridis_d(name = "Dataset", option = "turbo")
}

model <- function(k0, k1, k2, nsub, length) {
  k0 * sqrt(1 / (1 + k1 / nsub + k2 / (nsub * length)))
}

get_performance <- function(src = "data/performance") {
  d <- duckplyr::read_parquet_duckdb(
    fs::dir_ls(
      src,
      glob = "*parquet",
      recurse = TRUE
    ),
    prudence = "lavish"
  ) |>
    dplyr::summarise(
      r = mean(atanh(r)),
      .by = c(nsub, length, measure, iter, low)
    ) |>
    dplyr::summarise(r = tanh(mean(r)), .by = c(nsub, length, measure, low)) |>
    dplyr::mutate(
      length = length / 60,
      total = nsub * length,
      low = dplyr::replace_values(
        low,
        "False" ~ "High Motion",
        "True" ~ "Low Motion"
      )
    ) |>
    dplyr::collect()
}

test_performance <- function(performance) {
  d <- performance |>
    dplyr::filter(length == 58) |>
    dplyr::filter(total == max(total))

  lme4::lmer(atanh(r) ~ low + (1 | measure), data = d) |>
    emmeans::emmeans(~low) |>
    emmeans::contrast(method = "pairwise") |>
    broom::tidy(conf.int = TRUE)
}


get_key_measures <- function() {
  c(
    # these are the 18 that showed a good fit with their model
    "CardSort_Unadj",
    "Flanker_Unadj",
    "PMAT24_A_CR",
    "ReadEng_Unadj",
    "PicVocab_Unadj",
    "ProcSpeed_Unadj",
    "DDisc_AUC_40K",
    "VSPLOT_TC",
    "SCPT_SPEC",
    "PSQI_Score",
    "Emotion_Task_Face_Acc",
    "Language_Task_Story_Avg_Difficulty_Level",
    "Relational_Task_Acc",
    "WM_Task_Acc",
    "NEOFAC_O",
    "NEOFAC_C",
    "NEOFAC_E",
    "LifeSatisf_Unadj"
  )
}
