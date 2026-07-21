deg_2_rad <- function(rad) {
  rad * pi / 180
}

get_demographics <- function(datasets, by_run) {
  dplyr::bind_rows(datasets, .id = "dataset") |>
    dplyr::select(
      dataset,
      sub,
      sex,
      age,
      ses,
      bmi,
      deviceserialnumber,
      manufacturer
    ) |>
    dplyr::mutate(dataset = factor(dataset), sex = factor(sex)) |>
    dplyr::mutate(ses = dplyr::if_else(is.na(ses), "1", ses)) |>
    dplyr::distinct() |> # ABCD has some duplicates
    dplyr::mutate(
      sex = dplyr::replace_values(sex, "Male" ~ "M", "Female" ~ "F"),
      sub = stringr::str_remove(sub, "sub-"),
      ses = dplyr::if_else(is.na(ses), "1", ses)
    ) |>
    factorize_ses() |>
    dplyr::semi_join(
      dplyr::distinct(by_run, dataset, sub, ses),
      by = dplyr::join_by(dataset, sub, ses)
    )
}

factorize_ses <- function(.d) {
  .d |>
    dplyr::mutate(
      ses = factor(
        ses,
        levels = c("1", "2", "3", "4", "Baseline", "Year2", "Year4", "Year6"),
        ordered = TRUE
      )
    )
}

read_nda <- function(src) {
  header <- readr::read_tsv(here::here(src), n_max = 1, show_col_types = FALSE)

  readr::read_tsv(here::here(src), skip = 2, col_names = colnames(header))
}

write_demographics <- function(demographics, dst) {
  readr::write_tsv(demographics, dst)
  dst
}

summarise_by_time <- function(.data) {
  .data |>
    dplyr::summarise(
      sem = mad(framewise_displacement) / sqrt(dplyr::n()),
      median = median(framewise_displacement),
      .by = c(t, task, ses, scan, time, dataset, filtered)
    )
}


add_run <- function(.data) {
  if (!"run" %in% names(.data)) {
    .data <- .data |> dplyr::mutate(run = 1L)
  }
  .data
}

add_ped <- function(.data) {
  if (!"ped" %in% names(.data)) {
    .data <- .data |> dplyr::mutate(ped = NA_character_)
  }
  .data
}

do_casting <- function(.data) {
  .data |>
    dplyr::mutate(
      dplyr::across(tidyselect::any_of(c("sub", "ses")), as.character),
      dplyr::across(tidyselect::any_of(c("run", "scan")), as.integer)
    )
}


summarise_by <- function(.data, .cols) {
  .data |>
    dplyr::filter(t > 0) |>
    dplyr::summarise(
      sem = sd(framewise_displacement) / sqrt(dplyr::n()),
      loc = mean(framewise_displacement),
      .by = {{ .cols }}
    )
}

set_pedrun <- function(.data) {
  .data |>
    dplyr::mutate(
      ped2 = dplyr::replace_values(
        ped,
        "AP" ~ "2",
        "PA" ~ "1",
        "RL" ~ "1",
        "LR" ~ "2"
      ),
      ped.run = interaction(ped2, run)
    ) |>
    dplyr::select(-ped2)
}

.summarize_by_prop <- function(.data, threshold) {
  .data |>
    dplyr::summarise(
      lost = mean(framewise_displacement > threshold),
      .by = c(task, ses, sub, scan, dataset, filtered)
    )
}


.summarize_by_max <- function(.data, threshold) {
  .data |>
    dplyr::summarise(
      max_fd = max(framewise_displacement),
      .by = c(task, ses, sub, scan, dataset, filtered)
    )
}

get_lost_strict <- function(
  dataset,
  by_run,
  prop_thresh = 0.2,
  mfd_thresh = 0.25,
  max_fd_thresh = 5
) {
  n_subs <- by_run |>
    dplyr::semi_join(
      dplyr::distinct(dataset, dataset),
      by = dplyr::join_by(dataset)
    ) |>
    dplyr::count(dataset, task, ses, scan, filtered, name = "n_sub")

  by_prop <- .summarize_by_prop(dataset, threshold = 0.2) |>
    dplyr::filter(lost > prop_thresh) |>
    dplyr::distinct(dataset, task, ses, sub, scan, filtered)
  by_prop2 <- by_prop |>
    dplyr::count(dataset, task, ses, scan, filtered) |>
    dplyr::right_join(n_subs) |>
    dplyr::mutate(n = dplyr::if_else(is.na(n), 0, n)) |>
    dplyr::mutate(prop = n / n_sub) |>
    dplyr::select(-n, -n_sub)

  by_max <- .summarize_by_max(dataset) |>
    dplyr::filter(max_fd > max_fd_thresh) |>
    dplyr::distinct(dataset, task, ses, sub, scan, filtered)
  by_max2 <- by_max |>
    dplyr::count(dataset, task, ses, scan, filtered) |>
    dplyr::right_join(n_subs) |>
    dplyr::mutate(n = dplyr::if_else(is.na(n), 0, n)) |>
    dplyr::mutate(max = n / n_sub) |>
    dplyr::select(-n, -n_sub)

  by_avg <- by_run |>
    dplyr::semi_join(
      dplyr::distinct(dataset, dataset),
      by = dplyr::join_by(dataset)
    ) |>
    dplyr::filter(loc > mfd_thresh) |>
    dplyr::distinct(dataset, task, ses, sub, scan, filtered)
  by_avg2 <- by_avg |>
    dplyr::count(dataset, task, ses, scan, filtered) |>
    dplyr::right_join(n_subs) |>
    dplyr::mutate(n = dplyr::if_else(is.na(n), 0, n)) |>
    dplyr::mutate(avg = n / n_sub) |>
    dplyr::select(-n, -n_sub)

  dplyr::bind_rows(by_max, by_avg, by_prop) |>
    dplyr::distinct() |>
    dplyr::count(dataset, task, ses, scan, filtered) |>
    dplyr::right_join(
      n_subs,
      by = dplyr::join_by(dataset, task, ses, scan, filtered)
    ) |>
    dplyr::mutate(n = dplyr::if_else(is.na(n), 0, n)) |>
    dplyr::mutate(
      lost = n / n_sub,
      lower = qbeta(0.025, 1 / 2 + n, n_sub - n + 1 / 2),
      upper = qbeta(0.975, 1 / 2 + n, n_sub - n + 1 / 2)
    ) |>
    dplyr::left_join(
      by_prop2,
      by = dplyr::join_by(dataset, task, ses, scan, filtered)
    ) |>
    dplyr::left_join(
      by_max2,
      by = dplyr::join_by(dataset, task, ses, scan, filtered)
    ) |>
    dplyr::left_join(
      by_avg2,
      by = dplyr::join_by(dataset, task, ses, scan, filtered)
    )
}


get_lost_lenient <- function(by_run, mfd_thresh = 0.55) {
  n_subs <- by_run |>
    dplyr::count(dataset, task, ses, scan, filtered, name = "n_sub")

  by_run |>
    dplyr::filter(loc > mfd_thresh) |>
    dplyr::distinct(dataset, task, ses, sub, scan, filtered) |>
    dplyr::count(dataset, task, ses, scan, filtered) |>
    dplyr::right_join(n_subs) |>
    dplyr::mutate(n = dplyr::if_else(is.na(n), 0, n)) |>
    dplyr::mutate(
      lost = n / n_sub,
      lower = qbeta(0.025, 1 / 2 + n, n_sub - n + 1 / 2),
      upper = qbeta(0.975, 1 / 2 + n, n_sub - n + 1 / 2)
    )
}


bind_datasets <- function(datasets) {
  dplyr::bind_rows(datasets, .id = "dataset") |>
    factorize_ses() |>
    dplyr::mutate(scan = factor(scan, ordered = TRUE)) |>
    dplyr::select(
      dataset,
      sub,
      ses,
      task,
      scan,
      t,
      time,
      tidyselect::contains("frame")
    ) |>
    tidyr::pivot_longer(
      tidyselect::contains("frame"),
      names_to = "filtered",
      values_to = "framewise_displacement"
    ) |>
    dplyr::mutate(filtered = stringr::str_detect(filtered, "filtered")) |>
    dplyr::group_by(dataset)
}

bind_lost <- function(strict, lenient) {
  dplyr::bind_rows(list(lenient = lenient, strict = strict), .id = "type") |>
    dplyr::select(-n, -n_sub)
}


GeomSplitViolin <- ggplot2::ggproto(
  "GeomSplitViolin",
  ggplot2::GeomViolin,
  draw_group = function(self, data, ..., draw_quantiles = NULL) {
    data <- transform(
      data,
      xminv = x - violinwidth * (x - xmin),
      xmaxv = x + violinwidth * (xmax - x)
    )
    grp <- data[1, "group"]
    newdata <- plyr::arrange(
      transform(data, x = if (grp %% 2 == 1) xminv else xmaxv),
      if (grp %% 2 == 1) y else -y
    )
    newdata <- rbind(
      newdata[1, ],
      newdata,
      newdata[nrow(newdata), ],
      newdata[1, ]
    )
    newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[
      1,
      "x"
    ])
    if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
      stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <= 1))
      quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
      aesthetics <- data[
        rep(1, nrow(quantiles)),
        setdiff(names(data), c("x", "y")),
        drop = FALSE
      ]
      aesthetics$alpha <- rep(1, nrow(quantiles))
      both <- cbind(quantiles, aesthetics)
      quantile_grob <- GeomPath$draw_panel(both, ...)
      ggplot2:::ggname(
        "geom_split_violin",
        grid::grobTree(
          ggplot2::GeomPolygon$draw_panel(newdata, ...),
          quantile_grob
        )
      )
    } else {
      ggplot2:::ggname(
        "geom_split_violin",
        ggplot2::GeomPolygon$draw_panel(newdata, ...)
      )
    }
  }
)

geom_split_violin <- function(
  mapping = NULL,
  data = NULL,
  stat = "ydensity",
  position = "identity",
  ...,
  draw_quantiles = NULL,
  trim = TRUE,
  scale = "area",
  na.rm = FALSE,
  show.legend = NA,
  inherit.aes = TRUE
) {
  ggplot2::layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomSplitViolin,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      trim = trim,
      scale = scale,
      draw_quantiles = draw_quantiles,
      na.rm = na.rm,
      ...
    )
  )
}

write_png <- function(p, file, width, height) {
  ggplot2::ggsave(
    file,
    p,
    device = ragg::agg_png,
    width = width,
    height = height,
    create.dir = TRUE
  )
  file
}

write_tikz <- function(p, file, width, height) {
  ggplot2::ggsave(
    file,
    p,
    device = tikzDevice::tikz,
    width = width,
    height = height,
    standAlone = TRUE
  )
  file
}


get_lost_strict_for_modeling <- function(
  dataset,
  by_run,
  prop_thresh = 0.2,
  mfd_thresh = 0.25,
  max_fd_thresh = 5
) {
  by_prop <- .summarize_by_prop(dataset, threshold = 0.2) |>
    dplyr::mutate(exclude = lost > prop_thresh) |>
    dplyr::distinct(dataset, task, ses, sub, scan, filtered, exclude)

  by_max <- .summarize_by_max(dataset) |>
    dplyr::mutate(exclude = max_fd > max_fd_thresh) |>
    dplyr::distinct(dataset, task, ses, sub, scan, filtered, exclude)

  by_avg <- by_run |>
    dplyr::semi_join(
      dplyr::distinct(dataset, dataset),
      by = dplyr::join_by(dataset)
    ) |>
    dplyr::mutate(exclude = loc > mfd_thresh) |>
    dplyr::distinct(dataset, task, ses, sub, scan, filtered, exclude)

  dplyr::bind_rows(
    list(max = by_max, avg = by_avg, prop = by_prop),
    .id = "how"
  ) |>
    tidyr::pivot_wider(names_from = how, values_from = exclude) |>
    dplyr::mutate(exclude = max | avg | prop) |>
    dplyr::select(dataset, task, ses, sub, scan, filtered, exclude)
}

truncate_to_modal_lengths <- function(d) {
  expected <- d |>
    dplyr::summarise(n_tr = max(t), .by = c(sub, task, ses, scan)) |>
    dplyr::count(n_tr, task, ses, scan) |>
    dplyr::slice_max(
      order_by = n,
      n = 1,
      with_ties = FALSE,
      by = c(task, ses, scan)
    ) |>
    dplyr::select(-n)

  d |>
    dplyr::left_join(expected, by = dplyr::join_by(task, ses, scan)) |>
    dplyr::filter(t <= n_tr) |>
    dplyr::select(-n_tr)
}


# One combined demographics table across all datasets/sessions.
# Returns a tinytable, which renders natively to LaTeX/HTML/Typst/Word so the
# manuscript is not tied to any one output format. `format_tt(escape = TRUE)`
# is essential: the "%" in the Female column is a LaTeX comment character and
# will silently swallow row terminators if left unescaped.
make_demographics_table <- function(by_run, demographics) {
  fmt_mi <- function(x, digits = 0) {
    x <- x[!is.na(x)]
    if (!length(x)) {
      return(NA_character_)
    }
    q <- round(stats::quantile(x, c(0.5, 0.25, 0.75)), digits)
    glue::glue("{q[1]} ({q[2]}, {q[3]})")
  }
  ds_levels <- c(
    hcpya = "HCPYA",
    ukb = "UKB",
    hcpd = "HCPD",
    hcpa = "HCPA",
    abcd = "ABCD",
    spacetop = "SpaceTop"
  )
  demo <- by_run |>
    dplyr::distinct(sub, ses, dataset) |>
    dplyr::left_join(demographics, by = c("dataset", "sub", "ses")) |>
    dplyr::summarise(
      N = dplyr::n_distinct(sub),
      Age = fmt_mi(age, 1),
      BMI = fmt_mi(bmi, 1),
      Female = {
        f <- sum(sex == "F", na.rm = TRUE)
        sprintf("%s (%.0f%%)", format(f, big.mark = ","), 100 * f / dplyr::n())
      },
      .by = c(dataset, ses)
    ) |>
    dplyr::filter(!is.na(Age)) |>
    dplyr::mutate(Dataset = factor(ds_levels[dataset], levels = ds_levels)) |>
    dplyr::arrange(Dataset, ses) |>
    dplyr::mutate(
      N = format(N, big.mark = ","),
      BMI = tidyr::replace_na(BMI, "—"),
      Session = dplyr::if_else(ses == "1", "", as.character(ses)),
      Dataset = dplyr::if_else(duplicated(Dataset), "", as.character(Dataset))
    )

  demo |>
    dplyr::select(Dataset, Session, N, Age, BMI, Female) |>
    tinytable::tt() |>
    tinytable::format_tt(escape = TRUE) |>
    tinytable::style_tt(j = 1, bold = TRUE)
}

# High- vs low-mover demographics, one narrow tinytable. High/low are columns
# under per-dataset spanners so the two groups sit side by side; the group N is
# moved to its own row so the column headers can stay short ("Low"/"High") and
# the table fits the page width. See make_demographics_table for the escape note.
make_highlow_table <- function(group_highlow, demographics) {
  fmt_mi <- function(x, digits = 1) {
    x <- x[!is.na(x)]
    if (!length(x)) {
      return(NA_character_)
    }
    q <- round(stats::quantile(x, c(.5, .25, .75)), digits)
    sprintf("%s (%s, %s)", q[1], q[2], q[3])
  }
  pct <- function(sex, lvl) {
    n <- sum(sex == lvl, na.rm = TRUE)
    sprintf("%s (%.0f%%)", format(n, big.mark = ","), 100 * n / length(sex))
  }
  summ <- group_highlow |>
    dplyr::left_join(demographics, by = c("dataset", "sub")) |>
    dplyr::filter(
      (dataset == "abcd" & ses == "Baseline") | dataset == "hcpya"
    ) |>
    dplyr::summarise(
      N = format(dplyr::n_distinct(sub), big.mark = ","),
      Female = pct(sex, "F"),
      Male = pct(sex, "M"),
      Age = fmt_mi(age),
      BMI = fmt_mi(bmi),
      .by = c(dataset, split)
    ) |>
    dplyr::mutate(
      split = factor(split, levels = c("Low Movers", "High Movers")),
      dataset = factor(dataset, levels = c("abcd", "hcpya"))
    ) |>
    dplyr::arrange(dataset, split)

  wide <- summ |>
    dplyr::mutate(col = paste(dataset, split, sep = "|")) |>
    dplyr::select(col, N, Female, Male, Age, BMI) |>
    tidyr::pivot_longer(-col, names_to = "Characteristic") |>
    tidyr::pivot_wider(names_from = col, values_from = value)

  disp <- data.frame(
    Characteristic = wide$Characteristic,
    wide[["abcd|Low Movers"]],
    wide[["abcd|High Movers"]],
    wide[["hcpya|Low Movers"]],
    wide[["hcpya|High Movers"]],
    check.names = FALSE
  )
  names(disp) <- c("Characteristic", "Low", "High", "Low", "High")

  disp |>
    tinytable::tt() |>
    tinytable::format_tt(escape = TRUE) |>
    tinytable::group_tt(j = list("ABCD (Baseline)" = 2:3, "HCPYA" = 4:5))
}

read_distinct_collect <- function(source) {
  duckplyr::read_parquet_duckdb(source) |>
    dplyr::distinct(sub) |>
    dplyr::collect()
}

get_orig_counts <- function(
  ukb_source,
  spacetop_source,
  hcpya_source,
  hcpa_source,
  hcpd_source,
  abcd_source
) {
  abcd_entities <- get_abcd_entities(abcd_source) |>
    dplyr::distinct(sub) |>
    dplyr::mutate(dataset = "abcd")
  purrr::map(
    list(
      ukb = ukb_source,
      spacetop = spacetop_source,
      hcpya = hcpya_source,
      hcpa = hcpa_source,
      hcpd = hcpd_source
    ),
    read_distinct_collect
  ) |>
    dplyr::bind_rows(.id = "dataset") |>
    dplyr::bind_rows(abcd_entities) |>
    dplyr::count(dataset)
}

plot_spectrum <- function(spectrums) {
  ds <- unique(spectrums$dataset)
  d <- spectrums |>
    dplyr::mutate(avg = factor(avg, ordered = TRUE)) |>
    dplyr::select(-ses, -scan, -dataset, -sub)

  d |>
    ggplot2::ggplot(ggplot2::aes(y = avg, x = freq)) +
    ggplot2::facet_grid(task ~ param, scales = "free_y") +
    ggplot2::geom_raster(ggplot2::aes(fill = pxx)) +
    ggplot2::scale_fill_viridis_c(option = "turbo") +
    ggplot2::ylab("participant\n(<- lower avg fd, higher avg fd ->)") +
    ggplot2::theme(
      axis.ticks.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank()
    )
}

get_icc <- function(by_run) {
  ukb_abcd <- by_run |>
    dplyr::select(-sem) |>
    dplyr::filter(
      stringr::str_detect(task, "rest"),
      scan == 1,
      dataset %in% c("ukb", "abcd")
    ) |>
    dplyr::group_nest(dataset, filtered) |>
    dplyr::mutate(
      fit = purrr::map(
        data,
        ~ .x |>
          tidyr::pivot_wider(names_from = ses, values_from = loc) |>
          dplyr::select(-task, -scan, -sub) |>
          as.matrix() |>
          psych::ICC() |>
          purrr::pluck("results") |>
          tibble::as_tibble()
      )
    ) |>
    dplyr::select(-data) |>
    tidyr::unnest(fit)

  dd <- by_run |>
    dplyr::select(-sem) |>
    dplyr::filter(!(dataset %in% c("ukb", "spacetop"))) |>
    dplyr::filter(task == "rest") |>
    dplyr::group_nest(dataset, filtered, ses) |>
    dplyr::mutate(
      fit = purrr::map(
        data,
        ~ .x |>
          tidyr::pivot_wider(names_from = scan, values_from = loc) |>
          na.omit() |>
          dplyr::select(`1`:`4`) |>
          as.matrix() |>
          psych::ICC() |>
          purrr::pluck("results") |>
          tibble::as_tibble()
      )
    ) |>
    dplyr::select(-data) |>
    tidyr::unnest(fit) |>
    dplyr::bind_rows(ukb_abcd) |>
    dplyr::filter(type == "ICC1", !filtered) |>
    dplyr::rename(lower = `lower bound`, upper = `upper bound`)
}


format_icc <- function(row) {
  glue::glue(
    "$_{[round(row$lower, 2)]}{[round(row$ICC, 2)]}_{[round(row$upper, 2)]}$ ($F([row$df1],[row$df2]) = [round(row$F, 2)], p < 0.001$)",
    .open = "[",
    .close = "]"
  ) |>
    knitr::asis_output()
}


get_compare_datasets_fit <- function(by_run, demographics, bpm_src) {
  bpm <- duckplyr::read_parquet_duckdb(bpm_src) |>
    dplyr::select(
      sub = participant_id,
      ses = session_id,
      mh_t_bpm__ext_tscore,
      mh_t_bpm__int_tscore,
      mh_t_bpm__attn_tscore
    ) |>
    dplyr::filter(ses %in% c("ses-00A", "ses-02A", "ses-04A", "ses-06A")) |>
    dplyr::mutate(sub = stringr::str_remove(sub, "sub-")) |>
    convert_abcd_ses() |>
    tidyr::pivot_longer(mh_t_bpm__ext_tscore:mh_t_bpm__attn_tscore) |>
    na.omit() |>
    dplyr::summarise(value = mean(value), .by = c(sub, name)) |>
    tidyr::pivot_wider()

  d <- by_run |>
    dplyr::filter(
      dataset %in% c("hcpd", "abcd"),
      !filtered,
      stringr::str_detect(task, "rest")
    ) |>
    dplyr::left_join(dplyr::distinct(
      demographics,
      dataset,
      sub,
      age,
      ses,
      bmi,
      sex,
      deviceserialnumber,
      manufacturer
    )) |>
    dplyr::left_join(bpm) |>
    dplyr::mutate(dplyr::across(
      tidyselect::contains("mh_t_bpm"),
      ~ dplyr::if_else(is.na(.x), 0, .x)
    )) |>
    dplyr::mutate(
      external = mh_t_bpm__ext_tscore > 65,
      internal = mh_t_bpm__int_tscore > 65,
      attention = mh_t_bpm__attn_tscore > 65,
      deviceserialnumber = dplyr::if_else(
        is.na(deviceserialnumber),
        "hcpd",
        deviceserialnumber
      ),
      manufacturer = dplyr::if_else(
        is.na(manufacturer),
        "SIEMENS_HCD",
        manufacturer
      ),
      ses = factor(ses, ordered = FALSE)
    )

  glmmTMB::glmmTMB(
    loc ~ external +
      internal +
      attention +
      age * bmi * sex +
      manufacturer +
      manufacturer:age +
      (1 | sub),
    data = d
  )
}

format_tidy_contrast <- function(row, digits = 3) {
  if (row$p.value < 0.001) {
    p <- "p < 0.001"
  } else if (row$p.value < 0.01) {
    p <- "p < 0.01"
  } else if (row$p.value < 0.05) {
    p <- "p < 0.05"
  } else {
    p <- "p > 0.05"
  }
  str <- glue::glue(
    "$_{[round(row$conf.low, digits)]}{[round(row$estimate, digits)]}_{[round(row$conf.high, digits)]}$, $[p]$",
    .open = "[",
    .close = "]"
  )
  glue::glue("`[str]`{=latex}", .open = "[", .close = "]")
}


get_rest_fit <- function(by_run, demographics) {
  d <- by_run |>
    dplyr::filter(!filtered, stringr::str_detect(task, "rest")) |>
    dplyr::left_join(dplyr::distinct(
      demographics,
      dataset,
      sub,
      age,
      ses,
      bmi,
      sex,
    )) |>
    dplyr::mutate(ses = factor(ses, ordered = FALSE))

  glmmTMB::glmmTMB(
    loc ~ age * bmi * sex + dataset + dataset:age + (1 | sub),
    data = d
  )
}

plot_abcc_abcd <- function(abcd, mr_y_qc__mot) {
  abcc <- abcd |>
    dplyr::filter(task == "rest") |>
    dplyr::summarise(
      framewise_displacement = mean(framewise_displacement),
      .by = c(sub, ses)
    )

  abcd0 <- duckplyr::read_parquet_duckdb(mr_y_qc__mot, prudence = "lavish") |>
    dplyr::rename(sub = participant_id, ses = session_id) |>
    dplyr::mutate(sub = stringr::str_remove(sub, "sub-")) |>
    dplyr::select(
      sub,
      framewise_displacement = mr_y_qc__mot__rsfmri__mot_mean,
      ses
    ) |>
    dplyr::collect() |>
    convert_abcd_ses()

  dplyr::bind_rows(list(abcc = abcc, abcd = abcd0), .id = "dataset") |>
    tidyr::pivot_wider(
      names_from = dataset,
      values_from = framewise_displacement
    ) |>
    dplyr::mutate(.diff = abcd - abcc, .avg = (abcd + abcc) / 2) |>
    ggplot2::ggplot(ggplot2::aes(x = .avg, y = .diff)) +
    ggplot2::facet_wrap(~ses) +
    agree::geom_ba() +
    scattermore::geom_scattermore(pointsize = 2, alpha = 0.5) +
    ggplot2::ylab("Motion Difference: ABCD - ABCC (mm)") +
    ggplot2::xlab("Motion Average: (ABCD + ABCC)/2 (mm)")
}
