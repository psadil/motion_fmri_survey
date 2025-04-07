ledoit_wolf_shrinkage <- function(X){
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
  if (beta == 0){
    shrinkage <- 0
  }else{
    shrinkage <- beta / delta
  }
  shrinkage
}

lw <- function(X){
  shrinkage <- ledoit_wolf_shrinkage(X)
  n_samples <- nrow(X)
  n_features <- ncol(X)
  Xc <- scale(X, scale = FALSE)
  # emp_cov <- cov(X)
  emp_cov <- crossprod(Xc) / n_samples
  mu <- sum(psych::tr(emp_cov)) / n_features
  shrunk_cov <- (1.0 - shrinkage) * emp_cov
  shrunk_cov[seq(1, length(emp_cov), by=n_features+1)] <- shrunk_cov[seq(1, length(emp_cov), by=n_features+1)] + shrinkage * mu
  shrunk_cov
}


get_qc_src <- function(
    by_run, 
    timeseries_src="data/timeseries/derivatives", 
    quantile_range = c(.9, .99)
){
  avg <- by_run |>
    dplyr::filter(dataset=="ukb", ses=="2", task=="rest", !filtered) |>
    dplyr::distinct(sub, loc) |>
    dplyr::filter(
      dplyr::between(
        loc, 
        quantile(loc, quantile_range[1]),
        quantile(loc, quantile_range[2])
      )
    ) 
  
  
  tibble::tibble(src=fs::dir_ls(timeseries_src, glob="*arrow", recurse = TRUE)) |>
    dplyr::mutate(sub = stringr::str_extract(src, "[[:digit:]]{7}")) |>
    dplyr::semi_join(avg, by = dplyr::join_by(sub)) |>
    dplyr::left_join(avg) |>
    dplyr::arrange(dplyr::desc(loc)) |>
    dplyr::select(-sub, -loc)
}

get_qc_src_hcpya <- function(
    by_run,
    timeseries_src="/Users/psadil/Desktop/hcp_qc/timeseries.arrow", 
    quantile_range = c(.9, .99)
){
  avg <- by_run |>
    dplyr::filter(dataset=="hcpya", task=="rest", scan==2, !filtered) |>
    dplyr::distinct(sub, loc) |>
    dplyr::filter(
      dplyr::between(
        loc, 
        quantile(loc, quantile_range[1]),
        quantile(loc, quantile_range[2])
      )
    )
  
  tibble::tibble(src=fs::dir_ls(timeseries_src, glob="*arrow", recurse = TRUE)) |>
    dplyr::mutate(sub = stringr::str_extract(src, "[[:digit:]]{6}")) |>
    dplyr::semi_join(avg, by = dplyr::join_by(sub)) |>
    dplyr::left_join(avg) |>
    dplyr::arrange(dplyr::desc(loc)) |>
    dplyr::select(-sub, -loc)
}

get_cor_by_thresh <- function(d, window_width=150){
  
  n_tr <- dplyr::n_distinct(d$t)
  
  d |>
    tidyr::crossing(window_start = seq_len(n_tr - window_width) - 1) |>
    dplyr::select(-t) |>
    dplyr::group_nest(sub, ses, filtered, iter, window_start) |>
    dplyr::mutate(
      data = purrr::map2(
        data, window_start,
        ~ .x |>
          dplyr::slice_min(
            order_by = framewise_displacement,
            n=.y+window_width, 
            with_ties = FALSE
          ) |>
          dplyr::arrange(framewise_displacement) |>
          dplyr::slice_tail(n=window_width)
      ),
      max_fd = purrr::map_dbl(data, ~max(.x$framewise_displacement)),
      data = purrr::map(
        data, 
        ~ .x |>
          dplyr::select(-framewise_displacement) |>
          corrr::correlate(quiet = TRUE) |>
          corrr::shave() |> 
          corrr::stretch() |> 
          na.omit() 
      )
    ) |>
    tidyr::unnest(data)
}


get_cor_by_thresh_hcpya <- function(
    d, cleaned,
    gold,
    window_width=150, 
    timeseries_src="/Users/psadil/Desktop/hcp_qc/timeseries.arrow"){
  
  n_tr <- dplyr::n_distinct(d$t)
  
  ptseries <- arrow::open_dataset(timeseries_src, format = "ipc") |>
    dplyr::filter(sub==as.integer(unique(d$sub)), cleaned==.env$cleaned) |>
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
            n=.y+window_width, 
            with_ties = FALSE
          ) |>
          dplyr::arrange(framewise_displacement) |>
          dplyr::slice_tail(n=window_width)
      ),
      max_fd = purrr::map_dbl(data, ~max(.x$framewise_displacement)),
      data = purrr::map(
        data, 
        ~ .x |>
          dplyr::left_join(ptseries, by=dplyr::join_by(t)) |>
          dplyr::select(-t, -framewise_displacement) |>
          corrr::correlate(quiet = TRUE) |>
          corrr::shave() |> 
          corrr::stretch() |> 
          na.omit()
      )
    ) |>
    tidyr::unnest(data) |>
    dplyr::mutate(
      cleaned=.env$cleaned,
      x=stringr::str_extract(x, "[[:digit:]]+") |> as.integer(),
      y=stringr::str_extract(y, "[[:digit:]]+") |> as.integer()
    ) |>
    dplyr::left_join(gold, by = dplyr::join_by(x, y, filtered, sub, cleaned)) |>
    dplyr::summarise(
      r = cor(r.x, r.y),
      .by = c(window_start, max_fd, iter, filtered, sub, cleaned)
    )
}

get_cor_by_thresh_hcpya_gold <- function(
    d, cleaned,
    max_fd = 0.1,
    timeseries_src="/Users/psadil/Desktop/hcp_qc/timeseries.arrow"){
  
  d <- dplyr::filter(d, iter==0)
  n_tr <- dplyr::n_distinct(d$t)
  subs <- d |>
    dplyr::filter(iter==0) |>
    dplyr::distinct(sub) |>
    dplyr::mutate(sub = as.integer(sub))
  
  ptseries <- arrow::open_dataset(timeseries_src, format = "ipc") |>
    dplyr::filter(cleaned==.env$cleaned) |>
    dplyr::semi_join(subs, by = dplyr::join_by(sub)) |>
    dplyr::select(t, tidyselect::starts_with("column")) |>
    dplyr::collect()
  
  d_t <- d |>
    dplyr::filter(!filtered, iter==0) |>
    dplyr::filter(framewise_displacement < max_fd) |>
    dplyr::distinct(sub, t)
    
  d |>
    dplyr::semi_join(d_t, by = dplyr::join_by(sub, t)) |>
    dplyr::group_nest(sub, filtered, iter) |>
    dplyr::mutate(
      max_fd = purrr::map_dbl(data, ~max(.x$framewise_displacement)),
      mean_fd = purrr::map_dbl(data, ~mean(.x$framewise_displacement)),
      median_fd = purrr::map_dbl(data, ~median(.x$framewise_displacement)),
      data = purrr::map(
        data, 
        ~ .x |>
          dplyr::left_join(ptseries, by=dplyr::join_by(t)) |>
          dplyr::select(-t, -framewise_displacement) |>
          recipes::recipe() |>
          recipes::update_role(tidyselect::everything()) |>
          recipes::step_nzv(recipes::all_predictors()) |>
          recipes::prep() |>
          recipes::bake(new_data=NULL) |>
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
      cleaned=.env$cleaned,
      x=stringr::str_extract(x, "[[:digit:]]+") |> as.integer(),
      y=stringr::str_extract(y, "[[:digit:]]+") |> as.integer(),
      n_tr = nrow(d_t)
    )
}

get_qcs_gold_hcpya <- function(d){
  
  d |>
    dplyr::filter(window_start < 30, iter==0) |>
    dplyr::summarise(
      r = mean(r),
      .by = c(x, y, filtered, sub, cleaned)
    )
}

get_qcs_gold <- function(d){
  d |>
    dplyr::filter(window_start < 30, iter==0) |>
    dplyr::summarise(
      r = mean(r),
      .by = c(x, y, filtered, sub)
    )
}


get_cor_by_thresh2 <- function(qcs, gold){
  
  qcs |>
    dplyr::left_join(gold, by = dplyr::join_by(x, y, filtered, sub, ses)) |>
    dplyr::mutate(r = r.x - r.y) |>
    dplyr::summarise(
      r = median(r),
      .by = c(window_start, max_fd, iter, filtered, sub, ses)
    )
}

get_cor_by_thresh3 <- function(qcs, gold){
  
  qcs |>
    dplyr::left_join(gold, by = dplyr::join_by(x, y, filtered, sub, ses)) |>
    dplyr::summarise(
      r = cor(r.x, r.y),
      .by = c(window_start, max_fd, iter, filtered, sub, ses)
    )
}

get_cor_by_thresh3_hcpya <- function(qcs, gold){
  
  qcs |>
    dplyr::left_join(gold, by = dplyr::join_by(x, y, filtered, sub, cleaned)) |>
    dplyr::summarise(
      r = cor(r.x, r.y),
      .by = c(window_start, max_fd, iter, filtered, sub, cleaned)
    )
}


get_ukb_qc <- function(ukb){
  ukb |>
    dplyr::filter(ses=="2") |>
    dplyr::filter(stringr::str_detect(task, "rest")) |>
    dplyr::select(sub, ses, tidyselect::starts_with("frame"), t) |>
    tidyr::pivot_longer(
      tidyselect::starts_with("frame"), 
      names_to = "filtered", 
      values_to = "framewise_displacement") |>
    dplyr::mutate(filtered = stringr::str_detect(filtered, "filt", TRUE))
}

get_hcpya_qc <- function(hcpya){
  hcpya |>
    dplyr::filter(run==1) |>
    dplyr::filter(scan==1) |>
    dplyr::filter(stringr::str_detect(task, "rest")) |>
    dplyr::select(sub, tidyselect::starts_with("frame"), t) |>
    tidyr::pivot_longer(
      tidyselect::starts_with("frame"), 
      names_to = "filtered", 
      values_to = "framewise_displacement") |>
    dplyr::mutate(filtered = stringr::str_detect(filtered, "filt", TRUE))
}


get_qc_fd <- function(src, fd, shuffle, iter=0) {
  
  fd <- dplyr::filter(fd, sub == stringr::str_extract(src$src[[1]], "[[:digit:]]{7}"))
  
  d <- arrow::read_ipc_file(src$src[[1]]) |>
    dplyr::mutate(
      sub = stringr::str_extract(src, "(?<=sub=)[[:digit:]]+"),
      ses = stringr::str_extract(src, "(?<=ses=)[[:digit:]]+")
    ) |>
    dplyr::filter(t > 0) |>
    dplyr::left_join(fd, by = dplyr::join_by(t, sub, ses)) |>
    dplyr::select(-tar_group)
  
  if (shuffle){
    d <- d  |>
      dplyr::mutate(framewise_displacement = sample(framewise_displacement))
  }
  d |>
    dplyr::mutate(iter=iter, filtered=unique(fd$filtered))
}

get_qc_fd_hcpya <- function(hcpya, by_run, n_iter=100, n_sub=50, quantile_range = c(.5, .95)) {
  
  avg <- by_run |>
    dplyr::filter(dataset=="hcpya", task=="rest", scan==2, !filtered) |>
    dplyr::distinct(sub, loc) |>
    dplyr::filter(
      dplyr::between(
        loc, 
        quantile(loc, quantile_range[1]),
        quantile(loc, quantile_range[2])
      )
    ) |>
    dplyr::arrange(desc(loc)) |>
    dplyr::slice_max(order_by=loc, n=n_sub)
  
  fd <- hcpya |>
    dplyr::filter(task=="rest", scan==2) |>
    dplyr::semi_join(avg, by=dplyr::join_by(sub)) |>
    dplyr::select(sub, t, tidyselect::starts_with("frame")) |>
    tidyr::pivot_longer(
      tidyselect::starts_with("frame"), 
      names_to = "filtered",
      values_to = "framewise_displacement") |>
    dplyr::mutate(filtered = stringr::str_detect(filtered, "filt"))
    
  fd |> 
    tidyr::crossing(iter = seq_len(n_iter)) |>
    dplyr::group_nest(sub, iter, filtered) |>
    dplyr::mutate(
      data = purrr::map(
        data, 
        ~.x |>
          dplyr::mutate(framewise_displacement = sample(framewise_displacement)))
    ) |>
    tidyr::unnest(data) |>
    dplyr::bind_rows(dplyr::mutate(fd, iter=0))
}


bind_qc_fd <- function(...){
  dplyr::bind_rows(...)
}

get_qc_gold <- function(qcs){
  qcs |>
    dplyr::filter(window_start < 30, iter==0) |>
    dplyr::summarise(
      r = mean(r),
      .by = c(x, y, filtered, sub, ses)
    ) 
}

