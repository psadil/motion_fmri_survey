ledoit_wolf_shrinkage <- function(X) {
  n_samples <- nrow(X)
  n_features <- ncol(X)

  X <- scale(X, scale = FALSE)
  X2 <- X^2
  emp_cov_trace <- colSums(X2) / n_samples
  mu <- sum(emp_cov_trace) / n_features
  beta_ <- sum(crossprod(X2))
  delta_ <- sum(crossprod(X)^2) / n_samples^2
  # use delta_ to compute beta
  beta <- 1.0 / (n_features * n_samples) * (beta_ / n_samples - delta_)
  # delta is the sum of the squared coefficients of (<X.T,X> - mu*Id) / p
  delta <- (delta_ - 2.0 * mu * sum(emp_cov_trace) + n_features * mu^2) / n_features
  # get final beta as the min between beta and delta
  # We do this to prevent shrinking more than "1", which would invert
  # the value of covariances
  beta <- min(beta, delta)
  # finally get shrinkage
  if (beta == 0) {
    shrinkage <- 0
  } else {
    shrinkage <- beta / delta
  }
  shrinkage
}


lw <- function(X) {
  shrinkage <- ledoit_wolf_shrinkage(X)
  n_samples <- nrow(X)
  n_features <- ncol(X)
  Xc <- scale(X, scale = FALSE)
  # emp_cov <- cov(X)
  emp_cov <- crossprod(Xc) / n_samples
  mu <- sum(psych::tr(emp_cov)) / n_features
  shrunk_cov <- (1.0 - shrinkage) * emp_cov
  shrunk_cov[seq(1, length(emp_cov), by = n_features + 1)] <- shrunk_cov[seq(1, length(emp_cov), by = n_features + 1)] + shrinkage * mu
  shrunk_cov
}

get_cor_by_thresh_hcpya <- function(
    d, cleaned,
    gold,
    window_width = 50,
    timeseries_src = "data/timeseries.arrow") {
  n_tr <- dplyr::n_distinct(d$t)

  ptseries <- arrow::open_dataset(timeseries_src, format = "ipc") |>
    dplyr::filter(sub == as.integer(unique(d$sub)), cleaned == .env$cleaned) |>
    dplyr::select(t, tidyselect::starts_with("column")) |>
    dplyr::collect()

  d |>
    tidyr::crossing(window_start = seq_len(n_tr - window_width) - 1) |>
    dplyr::group_nest(sub, filtered, iter, window_start) |>
    dplyr::mutate(
      data = purrr::map2(
        data, window_start,
        ~ .x |>
          dplyr::slice_min(
            order_by = framewise_displacement,
            n = .y + window_width,
            with_ties = FALSE
          ) |>
          dplyr::arrange(framewise_displacement) |>
          dplyr::slice_tail(n = window_width)
      ),
      max_fd = purrr::map_dbl(data, ~ max(.x$framewise_displacement)),
      mean_fd = purrr::map_dbl(data, ~ mean(.x$framewise_displacement)),
      median_fd = purrr::map_dbl(data, ~ median(.x$framewise_displacement)),
      data = purrr::map(
        data,
        ~ .x |>
          dplyr::left_join(ptseries, by = dplyr::join_by(t)) |>
          dplyr::select(tidyselect::starts_with("column")) |>
          recipes::recipe() |>
          recipes::update_role(tidyselect::everything()) |>
          recipes::step_nzv(recipes::all_predictors()) |>
          recipes::prep() |>
          recipes::bake(new_data = NULL) |>
          as.matrix() |>
          lw() |>
          stats::cov2cor() |>
          corrr::as_cordf() |>
          # corrr::correlate(quiet = TRUE) |>
          corrr::shave() |>
          corrr::stretch() |>
          na.omit()
      )
    ) |>
    tidyr::unnest(data) |>
    dplyr::mutate(
      cleaned = .env$cleaned,
      x = stringr::str_extract(x, "[[:digit:]]+") |> as.integer(),
      y = stringr::str_extract(y, "[[:digit:]]+") |> as.integer()
    ) |>
    dplyr::left_join(gold, by = dplyr::join_by(x, y, filtered, sub, cleaned)) |>
    dplyr::summarise(
      r = cor(r.x, r.y),
      .by = c(window_start, iter, filtered, sub, cleaned, mean_fd, median_fd, max_fd)
    )
}

args <- commandArgs(trailingOnly = TRUE)
group <- as.integer(args[[1]])

qc_fd_hcpya_iter <- arrow::open_dataset("data/qc_fd_hcpya_iter.parquet") |>
  dplyr::filter(tar_group == group) |>
  dplyr::collect()
qcs_cor_hcpya_gold <- arrow::read_parquet("data/qcs_cor_hcpya_gold_max_0_lw.parquet") |>
  dplyr::select(sub, filtered, x, y, r, cleaned)

cleaned <- get_cor_by_thresh_hcpya(qc_fd_hcpya_iter, qcs_cor_hcpya_gold, cleaned = TRUE)
raw <- get_cor_by_thresh_hcpya(qc_fd_hcpya_iter, qcs_cor_hcpya_gold, cleaned = FALSE)

dplyr::bind_rows(cleaned, raw) |> arrow::write_parquet(glue::glue("derivatives/{group}.parquet"))
