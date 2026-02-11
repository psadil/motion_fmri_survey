library(dplyr)
library(ggplot2)
library(tidyr)
library(ggdist)
library(patchwork)
library(purrr)
library(gtsummary)

source("R/utils.R")


targets::tar_load(c(by_run, demographics))

mriqc <- arrow::open_dataset("data/bold") |>
  select(id = `_id`, fd_mean) |>
  collect() |>
  mutate(dataset = "mriqc")

all_by_run <- by_run |>
  filter(filtered) |>
  left_join(distinct(demographics, dataset, sub, ses, age, sex_gender, bmi)) |>
  mutate(
    dataset = replace_values(
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
  rename(sex = sex_gender)

ranges <- all_by_run |>
  summarise(
    xmin = quantile(age, 0.975, na.rm = TRUE),
    xmax = quantile(age, 0.025, na.rm = TRUE),
    .by = dataset
  ) |>
  mutate(y = seq(-.1, 0, by = 0.02))


a <- all_by_run |>
  ggplot(aes(x = age, color = dataset)) +
  scattermore::geom_scattermore(aes(y = loc), alpha = 0.5, pointsize = 1) +
  geom_segment(aes(x = xmin, xend = xmax, y = y), data = ranges) +
  geom_smooth(aes(y = loc, group = dataset), method = "lm") +
  coord_cartesian(ylim = c(-0.1, 1)) +
  scale_color_viridis_d(name = NULL, option = "turbo") +
  ylab("Framewise\nDisp. (mm)") +
  xlab("Age (years)") +
  theme_gray(base_size = 28) +
  guides(color = guide_legend(nrow = 1))

ranges_bmi <- all_by_run |>
  filter(stringr::str_detect(dataset, "SpaceTop", TRUE)) |>
  summarise(
    xmin = quantile(bmi, 0.975, na.rm = TRUE),
    xmax = quantile(bmi, 0.025, na.rm = TRUE),
    .by = dataset
  ) |>
  left_join(select(ranges, dataset, y), by = join_by(dataset))


b <- all_by_run |>
  filter(!is.na(bmi)) |>
  ggplot(aes(x = bmi, color = dataset)) +
  scattermore::geom_scattermore(aes(y = loc), alpha = 0.5, pointsize = 1) +
  geom_segment(aes(x = xmin, xend = xmax, y = y), data = ranges_bmi) +
  geom_smooth(aes(y = loc, group = dataset), method = "lm") +
  coord_cartesian(ylim = c(-0.1, 1)) +
  scale_color_viridis_d(option = "turbo", guide = "none", drop = FALSE) +
  ylab("Framewise\nDisp. (mm)") +
  xlab("Body Mass Index") +
  theme_gray(base_size = 28)

summary_data <- all_by_run |>
  filter(sex %in% c("F", "M")) |>
  summarise(
    mean = median(loc),
    min = quantile(loc, 0.25),
    max = quantile(loc, 0.75),
    .by = c(sex, dataset)
  )

c <- all_by_run |>
  filter(sex %in% c("F", "M")) |>
  ggplot(aes(x = dataset, color = sex, y = loc)) +
  geom_split_violin(trim = FALSE, alpha = 0.5) +
  geom_pointrange(
    data = summary_data,
    aes(dataset, mean, ymin = min, ymax = max, color = sex),
    shape = 20, # 95,
    position = position_dodge(width = 0.25)
  ) +
  coord_cartesian(ylim = c(-0.1, 1)) +
  ylab("Framewise\nDisp. (mm)") +
  xlab(NULL) +
  theme_gray(base_size = 28) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

d <- mriqc |>
  filter(fd_mean < quantile(fd_mean, 0.99)) |>
  ggplot(aes(x = dataset, y = fd_mean)) +
  stat_slab(
    aes(thickness = after_stat(pdf * n)),
    scale = 0.7,
    orientation = "vertical"
  ) +
  stat_dotsinterval(
    side = "bottom",
    scale = 0.7,
    slab_linewidth = NA,
    orientation = "vertical",
    quantiles = 100
  ) +
  coord_cartesian(ylim = c(-0.1, 1)) +
  ylab("Framewise\nDisp. (mm)") +
  xlab(NULL) +
  theme_gray(base_size = 28)

p <- a +
  b +
  c +
  d +
  plot_layout(nrow = 2, guides = "collect") &
  theme(legend.position = "bottom", legend.box = "vertical")

ggsave(
  "figures/age-chart.png",
  p,
  ragg::agg_png,
  width = 11,
  height = 8
)


p <- a +
  b2 +
  b1 +
  c +
  d +
  plot_layout(nrow = 1, guides = "collect", widths = c(3, 1, 3, 3, 1)) &
  theme(legend.position = "bottom")
ggsave(
  "figures/age-chart.png",
  p,
  width = 13,
  height = 6,
  device = ragg::agg_png
)


targets::tar_load(c(by_time, abcd_events, hcpya_events))

a2cps <- arrow::read_parquet("data/a2cpsbytime.parquet") |>
  mutate(time = t * .8)

by_time <- by_time |>
  dplyr::filter(!filtered)

plot_hcp_fd <- function(hcp_tmp, hcpya_events, task, fd_max_y = 0.5) {
  hcpya_events <- hcpya_events |>
    dplyr::mutate(scan = glue::glue("Scan: {scan}"))

  hcp_tmp <- hcp_tmp |>
    dplyr::mutate(scan = glue::glue("Scan: {scan}"))

  hcp_tmp |>
    filter(stringr::str_detect(task, .env$task)) |>
    ggplot(aes(x = time, y = median)) +
    facet_grid(scan ~ task, scales = "free_x") +
    geom_rect(
      data = filter(hcpya_events, stringr::str_detect(task, .env$task)),
      aes(
        xmin = onset,
        xmax = onset + duration,
        ymin = Inf,
        ymax = -Inf,
        fill = type
      ),
      alpha = 0.5,
      inherit.aes = FALSE
    ) +
    geom_line() +
    geom_ribbon(
      aes(ymin = median - 2 * sem, ymax = median + 2 * sem),
      alpha = 0.5
    ) +
    ylab("Framewise Disp. (mm)") +
    xlab("Time (s)") +
    scale_fill_viridis_d(option = "turbo") +
    theme_gray(base_size = 18) +
    theme(legend.position = "bottom") +
    coord_cartesian(ylim = c(0, fd_max_y))
}

hcp_tmp <- by_time |>
  filter(dataset == "hcpya")

a <- plot_hcp_fd(hcp_tmp, hcpya_events, "wm", fd_max_y = .2)
b <- plot_hcp_fd(hcp_tmp, hcpya_events, "gambling", fd_max_y = .2)
c <- plot_hcp_fd(hcp_tmp, hcpya_events, "motor", fd_max_y = .2)
d0 <- filter(by_time, dataset == "hcpya") |>
  filter(task == "rest") |>
  ggplot(aes(x = time, y = median, color = scan)) +
  facet_wrap(~task, scales = "free_x") +
  geom_line() +
  geom_ribbon(
    aes(ymin = median - 2 * sem, ymax = median + 2 * sem),
    alpha = 0.5
  ) +
  ylab("Framewise Disp. (mm)") +
  xlab("Time (s)") +
  theme_gray(base_size = 18) +
  theme(legend.position = "bottom") +
  coord_cartesian(ylim = c(0, 0.2))


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
  plot_layout(
    design = "
  1234
  5678
  "
  ) +
  plot_annotation("HCPYA (N: 1063)") &
  theme(legend.position = "bottom")

ggsave(
  "figures/hcp_by_time.png",
  p,
  ragg::agg_png,
  width = 13,
  height = 9
)

hcpd <- targets::tar_read(hcpd) |>
  count(task, run, ped, scan) |>
  slice_max(order_by = n, n = 1, by = c(task, scan)) |>
  select(task, ped, scan)

hcpd_events <- readr::read_csv(
  "data/events_consolidated/hcp_dev_event_times.csv"
) |>
  rename(ped = dir, type = trial_label) |>
  mutate(task = stringr::str_to_lower(task)) |>
  left_join(hcpd)

a <- plot_hcp_fd(
  filter(by_time, dataset == "hcpd") |> mutate(time = time / .72 * .8),
  hcpd_events,
  "carit",
  fd_max_y = .2
)
b <- plot_hcp_fd(
  filter(by_time, dataset == "hcpd") |> mutate(time = time / .72 * .8),
  hcpd_events,
  "emotion",
  fd_max_y = .2
)
c <- plot_hcp_fd(
  filter(by_time, dataset == "hcpd") |> mutate(time = time / .72 * .8),
  hcpd_events,
  "guessing",
  fd_max_y = .2
) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))

d <- filter(by_time, dataset == "hcpd") |>
  mutate(time = time / .72 * .8) |>
  filter(task == "rest") |>
  ggplot(aes(x = time, y = median, color = scan)) +
  facet_wrap(~task, scales = "free_x") +
  geom_line() +
  geom_ribbon(
    aes(ymin = median - 2 * sem, ymax = median + 2 * sem),
    alpha = 0.5
  ) +
  ylab("Framewise Disp. (mm)") +
  xlab("Time (s)") +
  guides(color = guide_legend(nrow = 2, byrow = TRUE)) +
  theme_gray(base_size = 20) +
  theme(legend.position = "bottom") +
  coord_cartesian(ylim = c(0, 0.2))

p <- a +
  b +
  c +
  d +
  plot_layout(nrow = 2) +
  plot_annotation(title = "HCPD (N: 652)") &
  theme(legend.position = "bottom")


ggsave(
  "figures/hcpd_by_time.png",
  p,
  ragg::agg_png,
  width = 6.5,
  height = 9
)


hcpa <- targets::tar_read(hcpa) |>
  count(task, run, ped, scan) |>
  slice_max(order_by = n, n = 1, by = c(task, scan)) |>
  select(task, ped, scan)

hcpa_events <- readr::read_csv(
  "data/events_consolidated/hcp_aging_event_times.csv"
) |>
  rename(ped = dir, type = trial_label) |>
  mutate(task = stringr::str_to_lower(task)) |>
  left_join(hcpa)

a <- plot_hcp_fd(
  filter(by_time, dataset == "hcpa") |>
    mutate(time = time / .72 * .8),
  hcpa_events,
  "carit",
  fd_max_y = .2
)
b <- plot_hcp_fd(
  filter(by_time, dataset == "hcpa") |> mutate(time = time / .72 * .8),
  hcpa_events,
  "facename",
  fd_max_y = .2
)
c <- plot_hcp_fd(
  filter(by_time, dataset == "hcpa") |> mutate(time = time / .72 * .8),
  hcpa_events,
  "vismotor",
  fd_max_y = .2
)

d <- filter(by_time, dataset == "hcpa") |>
  mutate(time = time / .72 * .8) |>
  filter(task == "rest") |>
  ggplot(aes(x = time, y = median, color = scan)) +
  facet_wrap(~task, scales = "free_x") +
  geom_line() +
  geom_ribbon(
    aes(ymin = median - 2 * sem, ymax = median + 2 * sem),
    alpha = 0.5
  ) +
  ylab("Framewise Disp. (mm)") +
  xlab("Time (s)") +
  guides(color = guide_legend(nrow = 2, byrow = TRUE)) +
  theme_gray(base_size = 20) +
  theme(legend.position = "bottom") +
  coord_cartesian(ylim = c(0, 0.2))


p <- a +
  b +
  c +
  d +
  plot_layout(nrow = 2) +
  plot_annotation(title = "HCPA (N: 725)") &
  theme(legend.position = "bottom")


ggsave(
  "figures/hcpa_by_time.png",
  p,
  ragg::agg_png,
  width = 6.5,
  height = 9
)

ukb_responses <- targets::tar_read(ukb_responses) |>
  mutate(task = "faces/shapes")


p <- by_time |>
  filter(dataset == "ukb", ses == "2") |>
  ggplot() +
  facet_wrap(~task, scales = "free_x", nrow = 1) +
  geom_vline(aes(xintercept = onset), data = ukb_responses, alpha = 0.5) +
  geom_rect(
    data = mutate(
      ukb_events,
      type = stringr::str_to_lower(type),
      task = "faces/shapes"
    ),
    aes(
      xmin = onset,
      xmax = onset + duration,
      fill = type
    ),
    ymax = Inf,
    ymin = -Inf,
    alpha = 0.2
  ) +
  geom_line(aes(x = time, y = median, color = scan)) +
  ylab("Framewise Disp. (mm)") +
  xlab("Time (s)") +
  coord_cartesian(ylim = c(0, 0.3)) +
  scale_colour_viridis_d(option = "turbo", drop = FALSE, guide = "none") +
  ggtitle("UKB (N: 42992)") +
  theme_gray(base_size = 20) +
  theme(legend.position = "bottom")

ggsave(
  "figures/ukb_by_time.png",
  p,
  ragg::agg_png,
  width = 13,
  height = 4
)

abcd_events <- readr::read_tsv("data/events_consolidated/abcd_events.tsv")

p <- by_time |>
  filter(dataset == "abcd") |>
  ggplot() +
  facet_grid(ses ~ task, scales = "free_x") +
  geom_line(aes(x = time, y = median, color = scan)) +
  ylab("Framewise Displacement") +
  xlab("Time (sec)") +
  coord_cartesian(ylim = c(0, 0.3)) +
  scale_colour_viridis_d(option = "turbo", drop = FALSE, guide = "none") +
  ggtitle("ABCD (N: 7929)") +
  geom_rect(
    aes(
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
  theme_gray(base_size = 20) +
  theme(legend.position = "bottom")

ggsave(
  "figures/abcd_by_time.png",
  p,
  ragg::agg_png,
  width = 13,
  height = 8
)


# spacetop_events <- readr::read_rds("data/events_consolidated/Spacetop_consolidated_events_table.rds") |>
#   as_tibble()

# sp_events <- readr::read_tsv("data/events/spacetop.tsv") |>
#
# spacetop_events +
#   rename(onset=median_onset_time, )
#   ggplot(aes(xmin=median_onset_time)) +
#   geom_rect(
#     aes(
#       xmin = onset,
#       xmax = onset + duration,
#       ymin = Inf,
#       ymax = -Inf,
#       fill = type
#     ),
#     alpha = 0.5,
#     inherit.aes = FALSE
#   )

p <- by_time |>
  filter(dataset == "spacetop", time < 400, time > 0) |>
  ggplot(aes(group = ses, color = scan)) +
  facet_wrap(~task, scales = "free") +
  geom_line(aes(x = time, y = median)) +
  ylab("Framewise Disp. (mm)") +
  xlab("Time (s)") +
  coord_cartesian(ylim = c(0, 0.2)) +
  ggtitle("SpaceTop (N: 113)") +
  scale_colour_viridis_d(option = "turbo", drop = FALSE) +
  guides(color = guide_legend(nrow = 1)) +
  theme_gray(base_size = 20) +
  theme(legend.position = "bottom")


ggsave(
  "figures/spacetop_by_time.png",
  p,
  ragg::agg_png,
  width = 13,
  height = 5
)


a2cps |>
  filter(t > 15) |>
  mutate(scan = factor(scan, ordered = TRUE)) |>
  ggplot() +
  facet_grid(~task, scales = "free_x") +
  geom_line(aes(x = time, y = framewise_displacement, color = scan)) +
  ylab("Framewise Displacement") +
  xlab("Time (sec)") +
  coord_cartesian(ylim = c(0, 0.3)) +
  scale_colour_viridis_d(option = "turbo", drop = FALSE, guide = "none") +
  ggtitle("A2CPS (N: 995)") +
  theme_gray(base_size = 20) +
  theme(legend.position = "bottom")


rest1_hcpya <- by_run |>
  filter(dataset == "hcpya", task == "rest", scan == 1, !filtered) |>
  mutate(split = loc > median(loc)) |>
  select(sub, split)

p <- by_run |>
  filter(dataset == "hcpya", !filtered) |>
  inner_join(rest1_hcpya) |>
  summarize(
    avg = mean(loc),
    sem = sd(loc) / sqrt(n()),
    .by = c(task, split, scan)
  ) |>
  mutate(split = if_else(split, "High Movers", "Low Movers")) |>
  ggplot(aes(x = task, y = avg, color = scan, group = scan)) +
  geom_point(
    aes(shape = split),
    position = position_dodge(width = 0.5),
    size = 3
  ) +
  geom_errorbar(
    aes(ymin = avg - 2 * sem, ymax = avg + 2 * sem),
    position = position_dodge(width = 0.5)
  ) +
  ylim(0, NA) +
  ylab("Framewise\nDisp. (mm)") +
  xlab(NULL) +
  theme_gray(base_size = 28) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1),
    legend.box = "vertical"
  )


ggsave(
  "figures/cluster-hcpya.png",
  p,
  ragg::agg_png,
  width = 10,
  height = 7
)


rest1_abcd <- by_run |>
  filter(
    dataset == "abcd",
    task == "rest",
    scan == 1,
    !filtered,
    ses == "baseline"
  ) |>
  mutate(split = loc > median(loc)) |>
  select(sub, split)

p <- by_run |>
  filter(dataset == "abcd", !filtered) |>
  inner_join(rest1_abcd) |>
  summarize(
    avg = mean(loc),
    sem = sd(loc) / sqrt(n()),
    .by = c(task, split, scan, ses)
  ) |>
  mutate(split = if_else(split, "High Movers", "Low Movers")) |>
  ggplot(aes(x = task, y = avg, color = scan, group = scan)) +
  facet_wrap(~ses, ncol = 1) +
  geom_point(
    aes(shape = split),
    position = position_dodge(width = 0.5),
    size = 3
  ) +
  geom_errorbar(
    aes(ymin = avg - 2 * sem, ymax = avg + 2 * sem),
    position = position_dodge(width = 0.5)
  ) +
  ylim(0, NA) +
  ylab("Framewise Disp. (mm)") +
  theme_gray(base_size = 28)


ggsave(
  "figures/cluster-abcd.png",
  p,
  ragg::agg_png,
  width = 10,
  height = 8
)


targets::tar_load(lost)

plot_lost <- function(lost, dataset, ses, base_size = 14) {
  p <- lost |>
    mutate(ses = glue::glue("ses: {ses}")) |>
    filter(!filtered, dataset == .env$dataset) |>
    ggplot(aes(y = task, fill = scan, x = lost)) +
    geom_col(aes(fill = scan), position = position_dodge(preserve = "single")) +
    geom_errorbarh(
      aes(xmin = lower, xmax = upper),
      position = position_dodge(preserve = "single"),
    ) +
    scale_x_continuous(NULL, breaks = c(0, 0.5), labels = c(0, 0.5)) +
    ylab(NULL) +
    guides(fill = guide_legend(nrow = 1)) +
    theme_gray(base_size = base_size)
  if (ses) {
    p <- p +
      facet_grid(ses ~ type)
  } else {
    p <- p +
      facet_wrap(~type, ncol = 2)
  }
  if (dataset == "spacetop") {
    p <- p +
      scale_fill_viridis_d(option = "turbo", drop = FALSE)
  } else {
    p <- p +
      scale_fill_viridis_d(option = "turbo", drop = FALSE) +
      guides(fill = "none")
  }
  p +
    ggtitle(dataset)
}

p <- plot_lost(lost, "abcd", TRUE) +
  plot_lost(lost, "hcpya", FALSE) +
  plot_lost(lost, "hcpd", FALSE) +
  plot_lost(lost, "hcpa", FALSE) +
  plot_lost(filter(lost, ses == 2), "ukb", FALSE) +
  plot_lost(lost, "spacetop", TRUE) +
  plot_layout(
    design = "
    23
    23
    23
    15
    16
    16
    46
    46
    46
    ",
    guides = "collect",
    ncol = 2
  ) &
  theme(legend.position = "bottom")


ggsave(
  "figures/all_motion_exclusion.png",
  p,
  device = ragg::agg_png,
  width = 6.5,
  height = 13
)

# lost by group

all_by_run |>
  filter(!is.na(age)) |>
  mutate(group = interaction(dataset, ses)) |>
  mutate(
    age = cut(age, breaks = quantile(age, c(0.25, 0.5, 0.75, 1))),
    N = n(),
    .by = group
  ) |>
  crossing(thresh = c(0.1, 0.2, 0.3, 0.4, 0.5, Inf)) |>
  filter(loc < thresh) |>
  mutate(thresh = factor(thresh, ordered = TRUE)) |>
  mutate(prop = n() / N, .by = c(group, thresh)) |>
  ggplot(aes(x = thresh, y = prop, color = age, group = age)) +
  geom_line() +
  facet_wrap(~group)

p <- all_by_run |>
  mutate(ses = forcats::fct_collapse(ses, " " = c("1", "2", "3", "4"))) |>
  group_nest(dataset, ses) |>
  mutate(
    data = purrr::map(
      data,
      ~ .x |>
        mutate(
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
  unnest(data) |>
  filter(!is.na(age), !is.na(sex)) |>
  mutate(N = n(), .by = c(dataset, ses, sex)) |>
  crossing(thresh = c(seq(0.1, 1, by = 0.1))) |>
  summarise(
    n = unique(sum(loc > thresh) / N),
    .by = c(dataset, thresh, age, ses, sex)
  ) |>
  ggplot(aes(x = thresh, y = n)) +
  geom_line(aes(color = age, linetype = sex)) +
  facet_wrap(~ dataset + ses, nrow = 1) +
  scale_color_viridis_d(option = "turbo", name = "Age Quartile") +
  xlab("Avg. Framewise Displacement Threshold") +
  ylab("Proportion Lost\nto Thresholding") +
  scale_x_continuous(
    "Avg. Framewise Displacement Threshold",
    limits = c(0, 1),
    breaks = c(0, 0.5),
    labels = c(0, 0.5)
  ) +
  theme_gray(base_size = 28) +
  theme(
    legend.position = "bottom"
  )

ggsave(
  "figures/motion-figs-2.png",
  p,
  device = ragg::agg_png,
  width = 14,
  height = 5
)

p <- all_by_run |>
  mutate(ses = forcats::fct_collapse(ses, " " = c("1", "2", "3", "4"))) |>
  filter(!is.na(bmi), !is.na(sex)) |>
  mutate(
    bmi = cut(
      bmi,
      breaks = c(0, 18.5, 25, 30, Inf),
      labels = c("Underweight", "Normal", "Overweight", "Obese"),
      ordered_result = TRUE
    )
  ) |>
  mutate(N = n(), .by = c(dataset, ses, sex)) |>
  crossing(thresh = c(seq(0.1, 1, by = 0.1))) |>
  summarise(
    n = unique(sum(loc > thresh) / N),
    .by = c(dataset, thresh, bmi, ses, sex)
  ) |>
  ggplot(aes(x = thresh, y = n)) +
  geom_line(aes(color = bmi, linetype = sex)) +
  facet_wrap(~ dataset + ses, nrow = 1) +
  scale_color_viridis_d(option = "turbo", name = "BMI") +
  ylab("Proportion Lost\nto Thresholding") +
  scale_x_continuous(
    "Avg. Framewise Displacement Threshold",
    limits = c(0, 1),
    breaks = c(0, 0.5),
    labels = c(0, 0.5)
  ) +
  theme_gray(base_size = 28) +
  theme(
    legend.position = "bottom"
  )

ggsave(
  "figures/motion-figs-2b.png",
  p,
  device = ragg::agg_png,
  width = 14,
  height = 5
)
