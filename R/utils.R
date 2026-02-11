deg_2_rad <- function(rad) {
  rad * pi / 180
}

get_demographics <- function(
  ukb,
  hcpya,
  spacetop,
  hcpdev,
  hcpaging,
  abcd
) {
  dplyr::bind_rows(
    list(
      ukb = ukb,
      hcpya = hcpya,
      spacetop = spacetop,
      hcpd = hcpdev,
      hcpa = hcpaging,
      abcd = abcd
    ),
    .id = "dataset"
  ) |>
    dplyr::mutate(
      scanner = dplyr::if_else(
        is.na(deviceserialnumber),
        site,
        deviceserialnumber
      )
    ) |>
    dplyr::select(
      dataset,
      sub,
      sex,
      gender,
      age,
      race,
      ethnicity,
      scanner,
      ses,
      bmi
    ) |>
    dplyr::mutate(
      dataset = factor(dataset),
      sex = factor(sex),
      gender = factor(gender),
      race = factor(race),
      ethnicity = factor(ethnicity),
      scanner = factor(scanner)
    ) |>
    dplyr::mutate(ses = dplyr::if_else(is.na(ses), "1", ses)) |>
    dplyr::distinct() |> # ABCD has some duplicates
    dplyr::mutate(
      sex_gender = dplyr::if_else(is.na(sex), gender, sex),
      sex_gender = dplyr::replace_values(
        sex_gender,
        "Male" ~ "M",
        "Female" ~ "F"
      ),
      sub = stringr::str_remove(sub, "sub-")
    ) |>
    dplyr::select(-sex, -gender) |>
    factorize_ses()
}

factorize_ses <- function(.d) {
  .d |>
    dplyr::mutate(
      ses = factor(
        ses,
        levels = c(
          "1",
          "2",
          "3",
          "4",
          "baseline",
          "Year2",
          "Year4"
        ),
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
    .data <- .data |>
      dplyr::mutate(run = 1L)
  }
  .data
}

add_ped <- function(.data) {
  if (!"ped" %in% names(.data)) {
    .data <- .data |>
      dplyr::mutate(ped = NA_character_)
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

get_by_run <- function(hcpya, hcpa, hcpd, ukb, abcd) {
  purrr::map(
    list(
      hcpya = dplyr::mutate(hcpya, ses = "1"),
      hcpa = hcpa,
      hcpd = hcpd,
      ukb = ukb,
      abcd = abcd
    ),
    summarise_by_run
  ) |>
    dplyr::bind_rows(.id = "dataset") |>
    dplyr::mutate(
      ses = factor(
        ses,
        levels = c(
          "1",
          "2",
          "3",
          "4",
          "baselineYear1Arm1",
          "2YearFollowUpYArm1",
          "4YearFollowUpYArm1"
        ),
        ordered = TRUE
      )
    ) |>
    set_pedrun()
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
    dplyr::semi_join(dplyr::distinct(dataset, dataset)) |>
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
    height = height
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
