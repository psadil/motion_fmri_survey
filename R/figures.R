make_fig_demographics <- function(
  by_run,
  demographics,
  mriqc_src,
  base_size = 11
) {
  mriqc <- arrow::open_dataset(mriqc_src) |>
    dplyr::select(id = `_id`, fd_mean) |>
    dplyr::collect() |>
    dplyr::mutate(dataset = "mriqc")

  all_by_run <- by_run |>
    dplyr::filter(filtered) |>
    dplyr::left_join(
      dplyr::distinct(
        demographics,
        dataset,
        sub,
        ses,
        age,
        sex_gender,
        bmi
      )
    ) |>
    dplyr::mutate(
      dataset = dplyr::replace_values(
        dataset,
        "abcd" ~ "ABCD",
        "hcpa" ~ "HCPA",
        "hcpd" ~ "HCPD",
        "hcpya" ~ "HCPYA",
        "ukb" ~ "UKB",
        "spacetop" ~ "SpaceTop"
      ) |>
        as.factor()
    ) |>
    dplyr::rename(sex = sex_gender)

  ranges <- all_by_run |>
    dplyr::summarise(
      xmin = quantile(age, 0.975, na.rm = TRUE),
      xmax = quantile(age, 0.025, na.rm = TRUE),
      .by = dataset
    ) |>
    dplyr::mutate(y = seq(-.1, 0, by = 0.02))

  a <- all_by_run |>
    ggplot2::ggplot(ggplot2::aes(x = age, color = dataset)) +
    scattermore::geom_scattermore(
      ggplot2::aes(y = loc),
      alpha = 0.5,
      pointsize = 1
    ) +
    ggplot2::geom_segment(
      ggplot2::aes(x = xmin, xend = xmax, y = y),
      data = ranges
    ) +
    ggplot2::geom_smooth(
      ggplot2::aes(y = loc, group = dataset),
      method = "lm"
    ) +
    ggplot2::coord_cartesian(ylim = c(-0.1, 1)) +
    ggplot2::scale_color_viridis_d(name = NULL, option = "turbo") +
    ggplot2::ylab("Framewise\nDisp. (mm)") +
    ggplot2::xlab("Age (years)") +
    ggplot2::theme_gray(base_size = base_size) +
    ggplot2::guides(color = ggplot2::guide_legend(nrow = 1))

  ranges_bmi <- all_by_run |>
    dplyr::filter(stringr::str_detect(dataset, "SpaceTop", TRUE)) |>
    dplyr::summarise(
      xmin = quantile(bmi, 0.975, na.rm = TRUE),
      xmax = quantile(bmi, 0.025, na.rm = TRUE),
      .by = dataset
    ) |>
    dplyr::left_join(
      dplyr::select(ranges, dataset, y),
      by = dplyr::join_by(dataset)
    )

  b <- all_by_run |>
    dplyr::filter(!is.na(bmi)) |>
    ggplot2::ggplot(ggplot2::aes(x = bmi, color = dataset)) +
    scattermore::geom_scattermore(
      ggplot2::aes(y = loc),
      alpha = 0.5,
      pointsize = 1
    ) +
    ggplot2::geom_segment(
      ggplot2::aes(x = xmin, xend = xmax, y = y),
      data = ranges_bmi
    ) +
    ggplot2::geom_smooth(
      ggplot2::aes(y = loc, group = dataset),
      method = "lm"
    ) +
    ggplot2::coord_cartesian(ylim = c(-0.1, 1)) +
    ggplot2::scale_color_viridis_d(
      option = "turbo",
      guide = "none",
      drop = FALSE
    ) +
    ggplot2::ylab("Framewise\nDisp. (mm)") +
    ggplot2::xlab("Body Mass Index") +
    ggplot2::theme_gray(base_size = base_size)

  summary_data <- all_by_run |>
    dplyr::filter(sex %in% c("F", "M")) |>
    dplyr::summarise(
      mean = median(loc),
      min = quantile(loc, 0.25),
      max = quantile(loc, 0.75),
      .by = c(sex, dataset)
    )

  c <- all_by_run |>
    dplyr::filter(sex %in% c("F", "M")) |>
    ggplot2::ggplot(ggplot2::aes(x = dataset, color = sex, y = loc)) +
    geom_split_violin(trim = FALSE, alpha = 0.5) +
    ggplot2::geom_pointrange(
      data = summary_data,
      ggplot2::aes(dataset, mean, ymin = min, ymax = max, color = sex),
      shape = 20, # 95,
      position = ggplot2::position_dodge(width = 0.25)
    ) +
    ggplot2::coord_cartesian(ylim = c(-0.1, 1)) +
    ggplot2::ylab("Framewise\nDisp. (mm)") +
    ggplot2::xlab(NULL) +
    ggplot2::theme_gray(base_size = base_size) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1)
    )

  d <- mriqc |>
    dplyr::filter(fd_mean < quantile(fd_mean, 0.99)) |>
    ggplot2::ggplot(ggplot2::aes(x = dataset, y = fd_mean)) +
    ggdist::stat_slab(
      ggplot2::aes(thickness = ggplot2::after_stat(pdf * n)),
      scale = 0.7,
      orientation = "vertical"
    ) +
    ggdist::stat_dotsinterval(
      side = "bottom",
      scale = 0.7,
      slab_linewidth = NA,
      orientation = "vertical",
      quantiles = 100
    ) +
    ggplot2::coord_cartesian(ylim = c(-0.1, 1)) +
    ggplot2::ylab("Framewise\nDisp. (mm)") +
    ggplot2::xlab(NULL) +
    ggplot2::theme_gray(base_size = base_size)

  a +
    b +
    c +
    d +
    patchwork::plot_layout(nrow = 2, guides = "collect") &
    ggplot2::theme(legend.position = "bottom", legend.box = "vertical")
}

plot_hcp_fd <- function(
  hcp_tmp,
  hcpya_events,
  task,
  fd_max_y = 0.5,
  base_size = 11
) {
  hcpya_events <- hcpya_events |>
    dplyr::mutate(scan = glue::glue("Scan: {scan}"))

  hcp_tmp <- hcp_tmp |>
    dplyr::mutate(scan = glue::glue("Scan: {scan}"))

  hcp_tmp |>
    dplyr::filter(stringr::str_detect(task, .env$task)) |>
    ggplot2::ggplot(ggplot2::aes(x = time, y = median)) +
    ggplot2::facet_grid(scan ~ task, scales = "free_x") +
    ggplot2::geom_rect(
      data = dplyr::filter(
        hcpya_events,
        stringr::str_detect(task, .env$task)
      ),
      ggplot2::aes(
        xmin = onset,
        xmax = onset + duration,
        ymin = Inf,
        ymax = -Inf,
        fill = type
      ),
      alpha = 0.5,
      inherit.aes = FALSE
    ) +
    ggplot2::geom_line() +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = median - 2 * sem, ymax = median + 2 * sem),
      alpha = 0.5
    ) +
    ggplot2::ylab("Framewise Disp. (mm)") +
    ggplot2::xlab("Time (s)") +
    ggplot2::scale_fill_viridis_d(option = "turbo") +
    ggplot2::theme_gray(base_size = base_size) +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::coord_cartesian(ylim = c(0, fd_max_y))
}

make_fig_by_time_hcpd <- function(
  hcpd,
  hcpd_event_times,
  by_time,
  filtered = FALSE,
  base_size = 11
) {
  hcpd <- hcpd |>
    dplyr::count(task, run, ped, scan) |>
    dplyr::slice_max(order_by = n, n = 1, by = c(task, scan)) |>
    dplyr::select(task, ped, scan)

  hcpd_events <- readr::read_csv(hcpd_event_times) |>
    dplyr::rename(ped = dir, type = trial_label) |>
    dplyr::mutate(task = stringr::str_to_lower(task)) |>
    dplyr::left_join(hcpd)

  a <- plot_hcp_fd(
    dplyr::filter(by_time, dataset == "hcpd", filtered == .env$filtered) |>
      dplyr::mutate(time = time / .72 * .8),
    hcpd_events,
    "carit",
    fd_max_y = .2,
    base_size = base_size
  )
  b <- plot_hcp_fd(
    dplyr::filter(by_time, dataset == "hcpd", filtered == .env$filtered) |>
      dplyr::mutate(time = time / .72 * .8),
    hcpd_events,
    "emotion",
    fd_max_y = .2,
    base_size = base_size
  )
  c <- plot_hcp_fd(
    dplyr::filter(by_time, dataset == "hcpd", filtered == .env$filtered) |>
      dplyr::mutate(time = time / .72 * .8),
    hcpd_events,
    "guessing",
    fd_max_y = .2,
    base_size = base_size
  ) +
    ggplot2::guides(fill = ggplot2::guide_legend(nrow = 1, byrow = TRUE))

  d <- dplyr::filter(by_time, dataset == "hcpd", filtered == .env$filtered) |>
    dplyr::mutate(time = time / .72 * .8) |>
    dplyr::filter(task == "rest") |>
    ggplot2::ggplot(ggplot2::aes(x = time, y = median, group = scan)) +
    ggplot2::facet_wrap(~task, scales = "free_x") +
    ggplot2::geom_line(ggplot2::aes(color = scan)) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = median - 2 * sem, ymax = median + 2 * sem),
      alpha = 0.5
    ) +
    ggplot2::ylab("Framewise Disp. (mm)") +
    ggplot2::xlab("Time (s)") +
    ggplot2::guides(color = ggplot2::guide_legend(nrow = 1, byrow = TRUE)) +
    ggplot2::theme_gray(base_size = base_size) +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::coord_cartesian(ylim = c(0, 0.2))

  a +
    b +
    c +
    d +
    patchwork::plot_layout(nrow = 2) +
    patchwork::plot_annotation(title = "HCPD")
}

make_fig_by_time_hcpa <- function(
  hcpa,
  hcpa_event_times,
  by_time,
  filtered = FALSE,
  base_size = 11
) {
  hcpa <- hcpa |>
    dplyr::count(task, run, ped, scan) |>
    dplyr::slice_max(order_by = n, n = 1, by = c(task, scan)) |>
    dplyr::select(task, ped, scan)

  hcpa_events <- readr::read_csv(hcpa_event_times, show_col_types = FALSE) |>
    dplyr::rename(ped = dir, type = trial_label) |>
    dplyr::mutate(task = stringr::str_to_lower(task)) |>
    dplyr::left_join(hcpa)

  a <- plot_hcp_fd(
    dplyr::filter(by_time, dataset == "hcpa", filtered == .env$filtered) |>
      dplyr::mutate(time = time / .72 * .8),
    hcpa_events,
    "carit",
    fd_max_y = .2,
    base_size = base_size
  )
  b <- plot_hcp_fd(
    dplyr::filter(by_time, dataset == "hcpa", filtered == .env$filtered) |>
      dplyr::mutate(time = time / .72 * .8),
    hcpa_events,
    "facename",
    fd_max_y = .2,
    base_size = base_size
  )
  c <- plot_hcp_fd(
    dplyr::filter(by_time, dataset == "hcpa", filtered == .env$filtered) |>
      dplyr::mutate(time = time / .72 * .8),
    hcpa_events,
    "vismotor",
    fd_max_y = 0.2,
    base_size = base_size
  )

  d <- dplyr::filter(by_time, dataset == "hcpa", filtered == .env$filtered) |>
    dplyr::mutate(time = time / .72 * .8) |>
    dplyr::filter(task == "rest") |>
    ggplot2::ggplot(ggplot2::aes(x = time, y = median, group = scan)) +
    ggplot2::facet_wrap(~task, scales = "free_x") +
    ggplot2::geom_line(ggplot2::aes(color = scan)) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = median - 2 * sem, ymax = median + 2 * sem),
      alpha = 0.5
    ) +
    ggplot2::ylab("Framewise Disp. (mm)") +
    ggplot2::xlab("Time (s)") +
    ggplot2::guides(color = ggplot2::guide_legend(nrow = 1, byrow = TRUE)) +
    ggplot2::theme_gray(base_size = base_size) +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::coord_cartesian(ylim = c(0, 0.2))

  p <- a +
    b +
    c +
    d +
    patchwork::plot_layout(nrow = 2) +
    patchwork::plot_annotation(title = "HCPA") &
    ggplot2::theme(legend.position = "bottom")
}

make_fig_by_time_hcpya <- function(
  by_time,
  hcpya_events,
  filtered = FALSE,
  base_size = 11
) {
  by_time <- by_time |>
    dplyr::filter(filtered == .env$filtered)

  hcp_tmp <- by_time |>
    dplyr::filter(dataset == "hcpya")

  a <- plot_hcp_fd(
    hcp_tmp,
    hcpya_events,
    "wm",
    fd_max_y = .2,
    base_size = base_size
  )
  b <- plot_hcp_fd(
    hcp_tmp,
    hcpya_events,
    "gambling",
    fd_max_y = .2,
    base_size = base_size
  )
  c <- plot_hcp_fd(
    hcp_tmp,
    hcpya_events,
    "motor",
    fd_max_y = .2,
    base_size = base_size
  )
  d0 <- dplyr::filter(by_time, dataset == "hcpya") |>
    dplyr::filter(task == "rest") |>
    ggplot2::ggplot(ggplot2::aes(x = time, y = median, group = scan)) +
    ggplot2::facet_wrap(~task, scales = "free_x") +
    ggplot2::geom_line(ggplot2::aes(color = scan)) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = median - 2 * sem, ymax = median + 2 * sem),
      alpha = 0.5
    ) +
    ggplot2::ylab("Framewise Disp. (mm)") +
    ggplot2::xlab("Time (s)") +
    ggplot2::theme_gray(base_size = base_size) +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::coord_cartesian(ylim = c(0, 0.2))

  d <- plot_hcp_fd(hcp_tmp, hcpya_events, "language", fd_max_y = .2)
  e <- plot_hcp_fd(hcp_tmp, hcpya_events, "social", fd_max_y = .2)
  f <- plot_hcp_fd(hcp_tmp, hcpya_events, "relational", fd_max_y = .2)
  g <- plot_hcp_fd(hcp_tmp, hcpya_events, "emotion", fd_max_y = .2)

  p <- a +
    b +
    c +
    d0 +
    d +
    e +
    f +
    g +
    patchwork::plot_layout(
      design = "
  1234
  5678
  "
    ) +
    patchwork::plot_annotation("HCPYA") &
    ggplot2::theme(legend.position = "bottom")
}

make_fig_by_time_ukb <- function(
  by_time,
  ukb_responses,
  ukb_events,
  filtered = FALSE,
  base_size = 11
) {
  by_time |>
    dplyr::filter(dataset == "ukb", ses == "2", filtered == .env$filtered) |>
    ggplot2::ggplot() +
    ggplot2::facet_wrap(~task, scales = "free_x", nrow = 1) +
    ggplot2::geom_vline(
      ggplot2::aes(xintercept = onset),
      data = ukb_responses,
      alpha = 0.5
    ) +
    ggplot2::geom_rect(
      data = dplyr::mutate(
        ukb_events,
        type = stringr::str_to_lower(type),
        task = "faces/shapes"
      ),
      ggplot2::aes(
        xmin = onset,
        xmax = onset + duration,
        fill = type
      ),
      ymax = Inf,
      ymin = -Inf,
      alpha = 0.2
    ) +
    ggplot2::geom_line(ggplot2::aes(x = time, y = median, color = scan)) +
    ggplot2::ylab("Framewise Disp. (mm)") +
    ggplot2::xlab("Time (s)") +
    ggplot2::coord_cartesian(ylim = c(0, 0.3)) +
    ggplot2::scale_colour_viridis_d(
      option = "turbo",
      drop = FALSE,
      guide = "none"
    ) +
    ggplot2::ggtitle("UKB") +
    ggplot2::theme_gray(base_size = base_size) +
    ggplot2::theme(legend.position = "bottom")
}

make_fig_by_time_abcd <- function(
  by_time,
  abcd_events_src,
  filtered = FALSE,
  base_size = 11
) {
  abcd_events <- readr::read_tsv(abcd_events_src, show_col_types = FALSE)
  by_time |>
    dplyr::filter(dataset == "abcd", filtered == .env$filtered) |>
    ggplot2::ggplot() +
    ggplot2::facet_grid(ses ~ task, scales = "free_x") +
    ggplot2::geom_line(ggplot2::aes(x = time, y = median, color = scan)) +
    ggplot2::ylab("Framewise Displacement") +
    ggplot2::xlab("Time (sec)") +
    ggplot2::coord_cartesian(ylim = c(0, 0.3)) +
    ggplot2::scale_colour_viridis_d(
      option = "turbo",
      drop = FALSE,
      guide = "none"
    ) +
    ggplot2::ggtitle("ABCD") +
    ggplot2::geom_rect(
      ggplot2::aes(
        xmin = onset,
        xmax = onset + duration,
        ymin = Inf,
        ymax = -Inf,
        fill = type
      ),
      data = abcd_events,
      alpha = 0.5,
      inherit.aes = FALSE
    ) +
    ggplot2::theme_gray(base_size = base_size) +
    ggplot2::theme(legend.position = "bottom")
}


make_fig_by_time_spacetop <- function(
  by_time,
  filtered = FALSE,
  base_size = 11
) {
  by_time |>
    dplyr::filter(
      dataset == "spacetop",
      time < 400,
      time > 0,
      filtered == .env$filtered
    ) |>
    ggplot2::ggplot(ggplot2::aes(group = ses, color = scan)) +
    ggplot2::facet_wrap(~task, scales = "free") +
    ggplot2::geom_line(ggplot2::aes(x = time, y = median)) +
    ggplot2::ylab("Framewise Disp. (mm)") +
    ggplot2::xlab("Time (s)") +
    ggplot2::coord_cartesian(ylim = c(0, 0.2)) +
    ggplot2::ggtitle("SpaceTop") +
    ggplot2::scale_colour_viridis_d(option = "turbo", drop = FALSE) +
    ggplot2::guides(color = ggplot2::guide_legend(nrow = 1)) +
    ggplot2::theme_gray(base_size = base_size) +
    ggplot2::theme(legend.position = "bottom")
}


make_fig_by_run <- function(by_run, filtered = FALSE, base_size = 10) {
  a <- by_run |>
    dplyr::filter(dataset == "hcpya", filtered == .env$filtered) |>
    ggplot2::ggplot(ggplot2::aes(y = task, color = scan, x = loc)) +
    ggplot2::geom_boxplot(
      outliers = FALSE,
      position = ggplot2::position_dodge(preserve = "single")
    ) +
    ggplot2::xlab("Framewise Disp. (mm)") +
    ggplot2::ylab(NULL) +
    ggplot2::coord_cartesian(xlim = c(0, 1)) +
    ggplot2::scale_colour_viridis_d(
      option = "turbo",
      drop = FALSE,
      guide = "none"
    ) +
    ggplot2::ggtitle("HCPYA")

  b <- by_run |>
    dplyr::filter(!is.na(scan), filtered == .env$filtered) |>
    dplyr::filter(dataset == "hcpd") |>
    ggplot2::ggplot(ggplot2::aes(y = task, color = scan, x = loc)) +
    ggplot2::geom_boxplot(
      outliers = FALSE,
      position = ggplot2::position_dodge(preserve = "single")
    ) +
    ggplot2::xlab("Framewise Disp. (mm)") +
    ggplot2::ylab(NULL) +
    ggplot2::coord_cartesian(xlim = c(0, 1)) +
    ggplot2::scale_colour_viridis_d(
      option = "turbo",
      drop = FALSE,
      guide = "none"
    ) +
    ggplot2::ggtitle("HCPD")

  c <- by_run |>
    dplyr::filter(dataset == "hcpa", filtered == .env$filtered) |>
    ggplot2::ggplot(ggplot2::aes(y = task, color = scan, x = loc)) +
    ggplot2::geom_boxplot(
      outliers = FALSE,
      position = ggplot2::position_dodge(preserve = "single")
    ) +
    ggplot2::xlab("Framewise Disp. (mm)") +
    ggplot2::ylab(NULL) +
    ggplot2::coord_cartesian(xlim = c(0, 1)) +
    ggplot2::scale_colour_viridis_d(
      option = "turbo",
      drop = FALSE,
      guide = "none"
    ) +
    ggplot2::ggtitle("HCPA")

  d <- by_run |>
    dplyr::filter(dataset == "ukb", filtered == .env$filtered) |>
    ggplot2::ggplot(ggplot2::aes(y = task, color = scan, x = loc)) +
    ggplot2::facet_wrap(~ses, nrow = 2, labeller = "label_both") +
    ggplot2::geom_boxplot(
      outliers = FALSE,
      position = ggplot2::position_dodge(preserve = "single")
    ) +
    ggplot2::xlab("Framewise Disp. (mm)") +
    ggplot2::ylab(NULL) +
    ggplot2::coord_cartesian(xlim = c(0, 1)) +
    ggplot2::scale_colour_viridis_d(
      option = "turbo",
      drop = FALSE,
      guide = "none"
    ) +
    ggplot2::ggtitle("UKB")

  e <- by_run |>
    dplyr::filter(!is.na(scan), filtered == .env$filtered) |>
    dplyr::filter(dataset == "abcd") |>
    ggplot2::ggplot(ggplot2::aes(y = task, color = scan, x = loc)) +
    ggplot2::facet_wrap(~ses, nrow = 3, labeller = "label_both") +
    ggplot2::geom_boxplot(
      outliers = FALSE,
      position = ggplot2::position_dodge(preserve = "single")
    ) +
    ggplot2::xlab("Framewise Disp. (mm)") +
    ggplot2::ylab(NULL) +
    ggplot2::coord_cartesian(xlim = c(0, 1)) +
    ggplot2::scale_colour_viridis_d(
      option = "turbo",
      drop = FALSE,
      guide = "none"
    ) +
    ggplot2::ggtitle("ABCD")

  f <- by_run |>
    dplyr::filter(!is.na(scan), filtered == .env$filtered) |>
    dplyr::filter(dataset == "spacetop") |>
    ggplot2::ggplot(ggplot2::aes(y = task, color = scan, x = loc)) +
    ggplot2::facet_wrap(~ses, scales = "free_y", labeller = "label_both") +
    ggplot2::geom_boxplot(
      outliers = FALSE,
      position = ggplot2::position_dodge(preserve = "single")
    ) +
    ggplot2::xlab("Framewise Disp. (mm)") +
    ggplot2::ylab(NULL) +
    ggplot2::coord_cartesian(xlim = c(0, 1)) +
    ggplot2::scale_colour_viridis_d(option = "turbo") +
    ggplot2::ggtitle("SpaceTop") +
    ggplot2::guides(color = ggplot2::guide_legend(nrow = 1))

  a +
    b +
    c +
    d +
    e +
    f +
    patchwork::plot_layout(
      guides = "collect",
      design = "
    555544
    555544
    223311
    666666
    666666
    "
    ) +
    patchwork::plot_annotation() &
    ggplot2::theme_gray(base_size = base_size) &
    ggplot2::theme(
      legend.position = "bottom"
    )
}

make_fig_all_motion_exclusion <- function(
  lost,
  filtered = FALSE,
  base_size = 7
) {
  a <- lost |>
    dplyr::filter(dataset == "hcpya", filtered == .env$filtered) |>
    ggplot2::ggplot(ggplot2::aes(y = task, x = lost)) +
    ggplot2::facet_wrap(~type) +
    ggplot2::geom_col(
      ggplot2::aes(fill = scan),
      position = ggplot2::position_dodge(preserve = "single")
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(xmin = lower, xmax = upper, group = scan),
      position = ggplot2::position_dodge(preserve = "single")
    ) +
    ggplot2::scale_fill_viridis_d(
      option = "turbo",
      drop = FALSE,
      guide = "none"
    ) +
    ggplot2::ggtitle("HCPYA") +
    ggplot2::theme_gray(base_size = base_size) +
    ggplot2::ylab(NULL) +
    ggplot2::xlab("Proportion Participants Excluded") +
    ggplot2::xlim(0, 1)

  b <- lost |>
    dplyr::filter(dataset == "hcpd", filtered == .env$filtered) |>
    ggplot2::ggplot(ggplot2::aes(y = task, x = lost)) +
    ggplot2::facet_wrap(~type) +
    ggplot2::geom_col(
      ggplot2::aes(fill = scan),
      position = ggplot2::position_dodge(preserve = "single")
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(xmin = lower, xmax = upper, group = scan),
      position = ggplot2::position_dodge(preserve = "single"),
    ) +
    ggplot2::scale_fill_viridis_d(
      option = "turbo",
      drop = FALSE,
      guide = "none"
    ) +
    ggplot2::ggtitle("HCPD") +
    ggplot2::theme_gray(base_size = base_size) +
    ggplot2::ylab(NULL) +
    ggplot2::xlab("Proportion Participants Excluded") +
    ggplot2::xlim(0, 1)

  c <- lost |>
    dplyr::filter(dataset == "hcpa", filtered == .env$filtered) |>
    ggplot2::ggplot(ggplot2::aes(y = task, x = lost)) +
    ggplot2::facet_wrap(~type) +
    ggplot2::geom_col(
      ggplot2::aes(fill = scan),
      position = ggplot2::position_dodge(preserve = "single")
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(xmin = lower, xmax = upper, group = scan),
      position = ggplot2::position_dodge(preserve = "single"),
    ) +
    ggplot2::scale_fill_viridis_d(
      option = "turbo",
      drop = FALSE,
      guide = "none"
    ) +
    ggplot2::ggtitle("HCPA") +
    ggplot2::theme_gray(base_size = base_size) +
    ggplot2::ylab(NULL) +
    ggplot2::xlab("Proportion Participants Excluded") +
    ggplot2::xlim(0, 1)

  d <- lost |>
    dplyr::filter(dataset == "ukb", filtered == .env$filtered) |>
    ggplot2::ggplot(ggplot2::aes(y = task, x = lost)) +
    ggplot2::facet_grid(ses ~ type) +
    ggplot2::geom_col(
      ggplot2::aes(fill = scan),
      position = ggplot2::position_dodge(preserve = "single")
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(xmin = lower, xmax = upper, group = scan),
      position = ggplot2::position_dodge(preserve = "single"),
    ) +
    ggplot2::scale_fill_viridis_d(
      option = "turbo",
      drop = FALSE,
      guide = "none"
    ) +
    ggplot2::ggtitle("UKB") +
    ggplot2::theme_gray(base_size = base_size) +
    ggplot2::ylab(NULL) +
    ggplot2::xlab("Proportion Participants Excluded") +
    ggplot2::xlim(0, 1)

  e <- lost |>
    dplyr::filter(dataset == "abcd", filtered == .env$filtered) |>
    ggplot2::ggplot(ggplot2::aes(y = task, x = lost)) +
    ggplot2::facet_grid(ses ~ type) +
    ggplot2::geom_col(
      ggplot2::aes(fill = scan),
      position = ggplot2::position_dodge(preserve = "single")
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(xmin = lower, xmax = upper, group = scan),
      position = ggplot2::position_dodge(preserve = "single"),
    ) +
    ggplot2::scale_fill_viridis_d(
      option = "turbo",
      drop = FALSE,
      guide = "none"
    ) +
    ggplot2::ggtitle("ABCD") +
    ggplot2::theme_gray(base_size = base_size) +
    ggplot2::ylab(NULL) +
    ggplot2::xlab("Proportion Participants Excluded") +
    ggplot2::xlim(0, 1)

  f <- lost |>
    dplyr::filter(dataset == "spacetop", filtered == .env$filtered) |>
    ggplot2::ggplot(ggplot2::aes(y = task, x = lost)) +
    ggplot2::facet_grid(ses ~ type, scales = "free_y") +
    ggplot2::geom_col(
      ggplot2::aes(fill = scan),
      position = ggplot2::position_dodge(preserve = "single")
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(xmin = lower, xmax = upper, group = scan),
      position = ggplot2::position_dodge(preserve = "single"),
    ) +
    ggplot2::scale_fill_viridis_d(option = "turbo", drop = FALSE) +
    ggplot2::ggtitle("Spacetop") +
    ggplot2::theme_gray(base_size = 7) +
    ggplot2::ylab(NULL) +
    ggplot2::xlab("Proportion Participants Excluded") +
    ggplot2::xlim(0, 1) +
    ggplot2::guides(fill = ggplot2::guide_legend(nrow = 1))

  a +
    b +
    c +
    d +
    e +
    f +
    patchwork::plot_layout(
      guides = "collect",
      design = "
    555111
    555111
    222333
    666444
    666444
    "
    ) &
    ggplot2::theme(legend.position = "bottom")
}

make_fig_cluster <- function(by_run, base_size = 10) {
  rest1_abcd <- by_run |>
    dplyr::filter(
      dataset == "abcd",
      task == "rest",
      scan == 1,
      !filtered,
      ses == "baseline"
    ) |>
    dplyr::mutate(split = loc > median(loc)) |>
    dplyr::select(sub, split)

  rest1_hcpya <- by_run |>
    dplyr::filter(dataset == "hcpya", task == "rest", scan == 1, !filtered) |>
    dplyr::mutate(split = loc > median(loc)) |>
    dplyr::select(sub, split)

  a <- by_run |>
    dplyr::filter(dataset == "hcpya", !filtered) |>
    dplyr::inner_join(rest1_hcpya) |>
    dplyr::summarize(
      avg = mean(loc),
      sem = sd(loc) / sqrt(dplyr::n()),
      .by = c(task, split, scan)
    ) |>
    dplyr::mutate(
      split = dplyr::if_else(split, "High Movers", "Low Movers") |>
        factor(levels = c("Low Movers", "High Movers"), ordered = TRUE)
    ) |>
    ggplot2::ggplot(
      ggplot2::aes(
        x = scan,
        y = avg,
        group = split
      )
    ) +
    ggplot2::facet_wrap(~task, scales = "free_x", nrow = 1) +
    ggplot2::geom_line(ggplot2::aes(linetype = split)) +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        ymin = avg - 2 * sem,
        ymax = avg + 2 * sem
      ),
      alpha = 0.5
    ) +
    ggplot2::ylim(0, NA) +
    ggplot2::ylab("Framewise\nDisp. (mm)") +
    ggplot2::xlab(NULL) +
    ggplot2::theme_gray(base_size = base_size) +
    ggplot2::ggtitle("HCPYA")

  b <- by_run |>
    dplyr::filter(dataset == "abcd", !filtered) |>
    dplyr::inner_join(rest1_abcd) |>
    dplyr::summarize(
      avg = mean(loc),
      sem = sd(loc) / sqrt(dplyr::n()),
      .by = c(task, split, scan, ses)
    ) |>
    dplyr::mutate(
      split = dplyr::if_else(split, "High Movers", "Low Movers") |>
        factor(levels = c("Low Movers", "High Movers"), ordered = TRUE)
    ) |>
    ggplot2::ggplot(
      ggplot2::aes(
        x = scan,
        y = avg,
        group = split
      )
    ) +
    ggplot2::facet_grid(ses ~ task, scales = "free_x") +
    ggplot2::geom_line(ggplot2::aes(linetype = split)) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = avg - 2 * sem, ymax = avg + 2 * sem),
      alpha = 0.5
    ) +
    ggplot2::ylim(0, NA) +
    ggplot2::ylab("Framewise Disp. (mm)") +
    ggplot2::theme_gray(base_size = base_size) +
    ggplot2::ggtitle("ABCD")

  a +
    b +
    patchwork::plot_layout(guides = "collect", ncol = 1, heights = c(1, 3)) +
    patchwork::plot_annotation(tag_levels = "a", tag_suffix = ")") &
    ggplot2::theme(
      legend.position = "bottom",
      legend.box = "vertical"
    )
}


make_fig_lost_by_group <- function(by_run, demographics, base_size = 10) {
  all_by_run <- by_run |>
    dplyr::filter(!filtered) |>
    dplyr::left_join(
      dplyr::distinct(
        demographics,
        dataset,
        sub,
        ses,
        age,
        sex_gender,
        bmi
      )
    ) |>
    dplyr::mutate(
      dataset = dplyr::replace_values(
        dataset,
        "abcd" ~ "ABCD",
        "hcpa" ~ "HCPA",
        "hcpd" ~ "HCPD",
        "hcpya" ~ "HCPYA",
        "ukb" ~ "UKB",
        "spacetop" ~ "SpaceTop"
      ) |>
        as.factor()
    ) |>
    dplyr::rename(sex = sex_gender)

  a <- all_by_run |>
    dplyr::mutate(
      ses = forcats::fct_collapse(ses, " " = c("1", "2", "3", "4"))
    ) |>
    dplyr::group_nest(dataset, ses) |>
    dplyr::mutate(
      data = purrr::map(
        data,
        ~ .x |>
          dplyr::mutate(
            age = cut(
              age,
              breaks = c(
                -Inf,
                quantile(age, probs = c(0.25, 0.5, 0.75, 1), na.rm = TRUE)
              ),
              labels = c("[0, 25)", "[25, 50)", "[50, 75)", "[75, 100]"),
              ordered_result = TRUE
            )
          )
      )
    ) |>
    tidyr::unnest(data) |>
    dplyr::filter(!is.na(age), !is.na(sex)) |>
    dplyr::mutate(N = dplyr::n(), .by = c(dataset, ses, sex)) |>
    tidyr::crossing(thresh = c(seq(0.1, 1, by = 0.1))) |>
    dplyr::summarise(
      n = unique(sum(loc > thresh) / N),
      .by = c(dataset, thresh, age, ses, sex)
    ) |>
    ggplot2::ggplot(ggplot2::aes(x = thresh, y = n)) +
    ggplot2::geom_line(ggplot2::aes(color = age, linetype = sex)) +
    ggplot2::facet_wrap(~ dataset + ses, nrow = 1) +
    ggplot2::scale_color_viridis_d(option = "turbo", name = "Age Quartile") +
    ggplot2::xlab("Avg. Framewise Displacement Threshold") +
    ggplot2::ylab("Proportion Lost\nto Thresholding") +
    ggplot2::scale_x_continuous(
      "Avg. Framewise Displacement Threshold",
      limits = c(0, 1),
      breaks = c(0, 0.5),
      labels = c(0, 0.5)
    ) +
    ggplot2::theme_gray(base_size = base_size) +
    ggplot2::theme(
      legend.position = "bottom"
    )

  b <- all_by_run |>
    dplyr::mutate(
      ses = forcats::fct_collapse(ses, " " = c("1", "2", "3", "4"))
    ) |>
    dplyr::filter(!is.na(bmi), !is.na(sex)) |>
    dplyr::mutate(
      bmi = cut(
        bmi,
        breaks = c(0, 18.5, 25, 30, Inf),
        labels = c("Underweight", "Normal", "Overweight", "Obese"),
        ordered_result = TRUE
      )
    ) |>
    dplyr::mutate(N = dplyr::n(), .by = c(dataset, ses, sex)) |>
    tidyr::crossing(thresh = c(seq(0.1, 1, by = 0.1))) |>
    dplyr::summarise(
      n = unique(sum(loc > thresh) / N),
      .by = c(dataset, thresh, bmi, ses, sex)
    ) |>
    ggplot2::ggplot(ggplot2::aes(x = thresh, y = n)) +
    ggplot2::geom_line(ggplot2::aes(color = bmi, linetype = sex)) +
    ggplot2::facet_wrap(~ dataset + ses, nrow = 1) +
    ggplot2::scale_color_viridis_d(option = "turbo", name = "BMI") +
    ggplot2::ylab("Proportion Lost\nto Thresholding") +
    ggplot2::scale_x_continuous(
      "Avg. Framewise Displacement Threshold",
      limits = c(0, 1),
      breaks = c(0, 0.5),
      labels = c(0, 0.5)
    ) +
    ggplot2::theme_gray(base_size = base_size) +
    ggplot2::theme(
      legend.position = "bottom"
    )

  a +
    b +
    patchwork::plot_layout(nrow = 2) +
    patchwork::plot_annotation(tag_levels = "a", tag_suffix = ")")
}
