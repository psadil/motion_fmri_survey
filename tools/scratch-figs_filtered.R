library(dplyr)
library(ggplot2)
library(tidyr)
library(ggdist)
library(patchwork)
library(purrr)
library(gtsummary)

targets::tar_load(c(by_run, demographics))
by_run <- by_run |>
  filter(filtered)

mriqc <- arrow::open_dataset("data/bold") |>
  select(id = `_id`, fd_mean) |>
  collect() |>
  mutate(dataset = "mriqc")

all_by_run <- by_run |>
  left_join(distinct(demographics, dataset, sub, ses, age, sex_gender, bmi)) |>
  mutate(
    dataset = case_match(
      dataset,
      "abcd" ~ "ABCD",
      "hcpa" ~ "HCPA",
      "hcpd" ~ "HCPD",
      "hcpya" ~ "HCPYA",
      "ukb" ~ "UKB",
      "spacetop" ~ "SpaceTop"
    ) |>
      as.factor()
  )

ranges <- all_by_run |>
  summarise(
    xmin = min(age, na.rm = TRUE),
    xmax = max(age, na.rm = TRUE),
    .by = dataset
  ) |>
  mutate(y = seq(-.1, 0, by = 0.02))

a <- all_by_run |>
  ggplot(aes(x = age, color = dataset)) +
  scattermore::geom_scattermore(aes(y = loc), alpha = 0.2, pointsize = 2) +
  geom_segment(aes(x = xmin, xend = xmax, y = y), data = ranges) +
  coord_cartesian(ylim = c(-0.1, 1)) +
  scale_color_viridis_d(option = "turbo") +
  ylab("Average Framewise\nDisplacement (mm)") +
  xlab("Age (years)") +
  theme_gray(base_size = 9) +
  guides(color = guide_legend(nrow = 1))

ranges_bmi <- all_by_run |>
  filter(stringr::str_detect(dataset, "ABCD|SpaceTop", TRUE)) |>
  summarise(
    xmin = min(bmi, na.rm = TRUE),
    xmax = max(bmi, na.rm = TRUE),
    .by = dataset
  ) |>
  left_join(select(ranges, dataset, y), by = join_by(dataset))


b <- all_by_run |>
  filter(!is.na(bmi)) |>
  ggplot(aes(x = bmi, color = dataset)) +
  scattermore::geom_scattermore(aes(y = loc), alpha = 0.2, pointsize = 2) +
  geom_segment(aes(x = xmin, xend = xmax, y = y), data = ranges_bmi) +
  coord_cartesian(ylim = c(-0.1, 1)) +
  scale_color_viridis_d(option = "turbo", guide = "none", drop = FALSE) +
  ylab("Average Framewise\nDisplacement (mm)") +
  xlab("Body Mass Index") +
  theme_gray(base_size = 9)

c <- all_by_run |>
  filter(sex_gender %in% c("F", "M")) |>
  ggplot(aes(x = sex_gender, color = dataset, y = loc)) +
  geom_boxplot(outliers = FALSE, position = position_dodge(preserve = "single")) +
  coord_cartesian(ylim = c(-0.1, 1)) +
  scale_color_viridis_d(option = "turbo", guide = "none", drop = FALSE) +
  ylab("Average Framewise\nDisplacement (mm)") +
  xlab("Sex/Gender") +
  theme_gray(base_size = 9)

d <- mriqc |>
  filter(fd_mean < quantile(fd_mean, 0.99)) |>
  ggplot(aes(x = dataset, y = fd_mean)) +
  stat_slab(aes(thickness = after_stat(pdf * n)), scale = 0.7, orientation = "vertical") +
  stat_dotsinterval(
    side = "bottom",
    scale = 0.7,
    slab_linewidth = NA,
    orientation = "vertical",
    quantiles = 100
  ) +
  coord_cartesian(ylim = c(-0.1, 1)) +
  ylab("Average Framewise\nDisplacement (mm)") +
  theme_gray(base_size = 9)

p <- a + b + c + d +
  plot_annotation(tag_levels = "a") +
  plot_layout(nrow = 2, guides = "collect") &
  theme(legend.position = "bottom")

ggsave(
  "figures/age-chart_filtered.png",
  p,
  width = 6.5,
  height = 4
)


a <- by_run |>
  filter(dataset == "hcpya") |>
  ggplot(aes(y = task, color = scan, x = loc)) +
  geom_boxplot(outliers = FALSE, position = position_dodge(preserve = "single")) +
  xlab("Framewise Displacement") +
  ylab(NULL) +
  coord_cartesian(xlim = c(0, 1)) +
  scale_colour_viridis_d(option = "turbo", drop = FALSE, guide = "none") +
  ggtitle("HCPYA")

b <- by_run |>
  filter(!is.na(scan)) |>
  filter(dataset == "hcpd") |>
  ggplot(aes(y = task, color = scan, x = loc)) +
  geom_boxplot(outliers = FALSE, position = position_dodge(preserve = "single")) +
  xlab("Framewise Displacement") +
  ylab(NULL) +
  coord_cartesian(xlim = c(0, 1)) +
  scale_colour_viridis_d(option = "turbo", drop = FALSE, guide = "none") +
  ggtitle("HCPD")

c <- by_run |>
  filter(dataset == "hcpa") |>
  ggplot(aes(y = task, color = scan, x = loc)) +
  geom_boxplot(outliers = FALSE, position = position_dodge(preserve = "single")) +
  xlab("Framewise Displacement") +
  ylab(NULL) +
  coord_cartesian(xlim = c(0, 1)) +
  scale_colour_viridis_d(option = "turbo", drop = FALSE, guide = "none") +
  ggtitle("HCPA")

d <- by_run |>
  filter(dataset == "ukb") |>
  ggplot(aes(y = task, color = scan, x = loc)) +
  facet_wrap(~ses, nrow = 2) +
  geom_boxplot(outliers = FALSE, position = position_dodge(preserve = "single")) +
  xlab("Framewise Displacement") +
  ylab(NULL) +
  coord_cartesian(xlim = c(0, 1)) +
  scale_colour_viridis_d(option = "turbo", drop = FALSE, guide = "none") +
  ggtitle("UKB")

e <- by_run |>
  filter(!is.na(scan)) |>
  mutate(scan = factor(scan)) |>
  filter(dataset == "abcd") |>
  mutate(
    ses = forcats::fct_recode(
      ses,
      baseline = "baselineYear1Arm1",
      Year2 = "2YearFollowUpYArm1",
      Year4 = "4YearFollowUpYArm1"
    )
  ) |>
  ggplot(aes(y = task, color = scan, x = loc)) +
  facet_wrap(~ses, nrow = 3) +
  geom_boxplot(outliers = FALSE, position = position_dodge(preserve = "single")) +
  xlab("Framewise Displacement") +
  ylab(NULL) +
  coord_cartesian(xlim = c(0, 1)) +
  scale_colour_viridis_d(option = "turbo", drop = FALSE, guide = "none") +
  ggtitle("ABCD")

f <- by_run |>
  filter(!is.na(scan)) |>
  mutate(scan = factor(scan)) |>
  filter(dataset == "spacetop") |>
  ggplot(aes(y = task, color = scan, x = loc)) +
  facet_wrap(~ses, scales = "free_y") +
  geom_boxplot(outliers = FALSE, position = position_dodge(preserve = "single")) +
  xlab("Framewise Displacement") +
  ylab(NULL) +
  coord_cartesian(xlim = c(0, 1)) +
  scale_colour_viridis_d(option = "turbo") +
  ggtitle("SpaceTop") +
  guides(color = guide_legend(nrow = 1))

p <- a + b + c + d + e + f +
  plot_layout(
    guides = "collect",
    design = "
    555544
    223311
    666666
    "
  ) +
  plot_annotation() &
  theme_gray(base_size = 7) &
  theme(
    legend.position = "bottom"
  )

ggsave(
  "figures/all_motion-avg_filtered.png",
  p,
  device = ragg::agg_png,
  width = 6.5,
  height = 7.5
)


targets::tar_load(c(by_time, abcd_events, hcpya_events))
by_time <- by_time |>
  filter(filtered)

hcpya_events <- hcpya_events |>
  mutate(
    task = forcats::fct_relabel(task, stringr::str_to_lower),
    ymin = if_else(ped == "LR", 0, 0.15),
    ymax = if_else(ped == "LR", 0.15, 0.3)
  ) |>
  filter(stringr::str_detect(task, "language", TRUE))

a <- by_time |>
  filter(dataset == "hcpya") |>
  ggplot() +
  facet_wrap(~task, scales = "free_x") +
  geom_line(aes(x = time, y = median, color = scan)) +
  ylab("Framewise Displacement") +
  xlab("Time (sec)") +
  coord_cartesian(ylim = c(0, 0.3)) +
  scale_colour_viridis_d(option = "turbo", drop = FALSE, guide = "none") +
  ggtitle("HCPYA") +
  geom_rect(
    aes(
      xmin = onset,
      xmax = onset + duration,
      ymin = ymin,
      ymax = ymax,
      fill = type
    ),
    data = hcpya_events,
    alpha = 0.5,
    inherit.aes = FALSE,
    show.legend = FALSE
  )

b <- by_time |>
  filter(dataset == "hcpd") |>
  ggplot() +
  facet_wrap(~task, scales = "free_x") +
  geom_line(aes(x = time, y = median, color = scan)) +
  ylab("Framewise Displacement") +
  xlab("Time (sec)") +
  coord_cartesian(ylim = c(0, 0.3)) +
  scale_colour_viridis_d(option = "turbo", drop = FALSE, guide = "none") +
  ggtitle("HCPD")

c <- by_time |>
  filter(dataset == "hcpa") |>
  ggplot() +
  facet_wrap(~task, scales = "free_x") +
  geom_line(aes(x = time, y = median, color = scan)) +
  ylab("Framewise Displacement") +
  xlab("Time (sec)") +
  coord_cartesian(ylim = c(0, 0.3)) +
  scale_colour_viridis_d(option = "turbo", drop = FALSE, guide = "none") +
  ggtitle("HCPA")

d <- by_time |>
  filter(dataset == "ukb", ses == "2") |>
  ggplot() +
  facet_wrap(~task, scales = "free_x", nrow = 2) +
  geom_line(aes(x = time, y = median, color = scan)) +
  ylab("Framewise Displacement") +
  xlab("Time (sec)") +
  coord_cartesian(ylim = c(0, 0.3)) +
  scale_colour_viridis_d(option = "turbo", drop = FALSE, guide = "none") +
  ggtitle("UKB")

e <- by_time |>
  filter(dataset == "abcd") |>
  ggplot() +
  facet_grid(ses ~ task, scales = "free_x") +
  geom_line(aes(x = time, y = median, color = scan)) +
  ylab("Framewise Displacement") +
  xlab("Time (sec)") +
  coord_cartesian(ylim = c(0, 0.3)) +
  scale_colour_viridis_d(option = "turbo", drop = FALSE, guide = "none") +
  ggtitle("ABCD") +
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
    inherit.aes = FALSE,
    show.legend = FALSE
  )

f <- by_time |>
  filter(dataset == "spacetop", time < 400) |>
  ggplot(aes(group = ses)) +
  facet_wrap(~task, scales = "free") +
  geom_line(aes(x = time, y = median, color = scan)) +
  ylab("Framewise Displacement") +
  xlab("Time (sec)") +
  coord_cartesian(ylim = c(0, 0.3)) +
  scale_colour_viridis_d(option = "turbo", drop = FALSE, guide = "none") +
  ggtitle("SpaceTop") +
  guides(color = guide_legend(nrow = 1))

p <- a + b + c + d + e + f +
  plot_layout(
    guides = "collect",
    design = "
    111111
    222333
    444555
    666666
    "
  ) +
  plot_annotation(tag_levels = "a") &
  theme_gray(base_size = 7) &
  theme(
    legend.position = "bottom"
  )


ggsave(
  "figures/all_motion_filtered.png",
  p,
  device = ragg::agg_png,
  width = 6.5,
  height = 9.5
)

ggsave(
  "figures/all_motion.tex",
  p,
  tikzDevice::tikz,
  width = 6.5,
  height = 7.5
)


lost <- targets::tar_read(lost) |>
  filter(filtered)


a <- lost |>
  filter(dataset == "hcpya") |>
  ggplot(aes(y = task, x = lost)) +
  facet_wrap(~type) +
  geom_col(aes(fill = scan), position = position_dodge(preserve = "single")) +
  geom_errorbarh(
    aes(xmin = lower, xmax = upper, group = scan),
    position = position_dodge(preserve = "single"),
  ) +
  scale_fill_viridis_d(option = "turbo", drop = FALSE, guide = "none") +
  ggtitle("HCPYA") +
  theme_gray(base_size = 7) +
  ylab(NULL) +
  xlab("Proportion Participants Excluded") +
  xlim(0, 1)

b <- lost |>
  filter(dataset == "hcpd") |>
  ggplot(aes(y = task, x = lost)) +
  facet_wrap(~type) +
  geom_col(aes(fill = scan), position = position_dodge(preserve = "single")) +
  geom_errorbarh(
    aes(xmin = lower, xmax = upper, group = scan),
    position = position_dodge(preserve = "single"),
  ) +
  scale_fill_viridis_d(option = "turbo", drop = FALSE, guide = "none") +
  ggtitle("HCPD") +
  theme_gray(base_size = 7) +
  ylab(NULL) +
  xlab("Proportion Participants Excluded") +
  xlim(0, 1)

c <- lost |>
  filter(dataset == "hcpa") |>
  ggplot(aes(y = task, x = lost)) +
  facet_wrap(~type) +
  geom_col(aes(fill = scan), position = position_dodge(preserve = "single")) +
  geom_errorbarh(
    aes(xmin = lower, xmax = upper, group = scan),
    position = position_dodge(preserve = "single"),
  ) +
  scale_fill_viridis_d(option = "turbo", drop = FALSE, guide = "none") +
  ggtitle("HCPA") +
  theme_gray(base_size = 7) +
  ylab(NULL) +
  xlab("Proportion Participants Excluded") +
  xlim(0, 1)


d <- lost |>
  filter(dataset == "ukb") |>
  ggplot(aes(y = task, x = lost)) +
  facet_grid(ses ~ type) +
  geom_col(aes(fill = scan), position = position_dodge(preserve = "single")) +
  geom_errorbarh(
    aes(xmin = lower, xmax = upper, group = scan),
    position = position_dodge(preserve = "single"),
  ) +
  scale_fill_viridis_d(option = "turbo", drop = FALSE, guide = "none") +
  ggtitle("UKB") +
  theme_gray(base_size = 7) +
  ylab(NULL) +
  xlab("Proportion Participants Excluded") +
  xlim(0, 1)

e <- lost |>
  filter(dataset == "abcd") |>
  ggplot(aes(y = task, x = lost)) +
  facet_grid(ses ~ type) +
  geom_col(aes(fill = scan), position = position_dodge(preserve = "single")) +
  geom_errorbarh(
    aes(xmin = lower, xmax = upper, group = scan),
    position = position_dodge(preserve = "single"),
  ) +
  scale_fill_viridis_d(option = "turbo", drop = FALSE, guide = "none") +
  ggtitle("ABCD") +
  theme_gray(base_size = 7) +
  ylab(NULL) +
  xlab("Proportion Participants Excluded") +
  xlim(0, 1)

f <- lost |>
  filter(dataset == "spacetop") |>
  ggplot(aes(y = task, x = lost)) +
  facet_grid(ses ~ type, scales = "free_y") +
  geom_col(aes(fill = scan), position = position_dodge(preserve = "single")) +
  geom_errorbarh(
    aes(xmin = lower, xmax = upper, group = scan),
    position = position_dodge(preserve = "single"),
  ) +
  scale_fill_viridis_d(option = "turbo", drop = FALSE) +
  ggtitle("Spacetop") +
  theme_gray(base_size = 7) +
  ylab(NULL) +
  xlab("Proportion Participants Excluded") +
  xlim(0, 1) +
  guides(fill = guide_legend(nrow = 1))


p <- a + b + c + d + e + f +
  plot_layout(
    guides = "collect",
    design = "
    555111
    555111
    222333
    666444
    666444
    "
  ) &
  theme(legend.position = "bottom")

ggsave(
  "figures/all_motion_exclusion_filtered.png",
  p,
  device = ragg::agg_png,
  width = 6.5,
  height = 6.5
)


lost2 <- targets::tar_read(lost_strict) |>
  filter(filtered) |>
  select(dataset:scan, prop, max, avg) |>
  pivot_longer(
    c(prop, max, avg),
    names_to = "type",
    values_to = "lost"
  )

a <- lost2 |>
  filter(dataset == "hcpya") |>
  ggplot(aes(y = task, x = lost)) +
  facet_wrap(~type) +
  geom_col(aes(fill = scan), position = position_dodge(preserve = "single")) +
  scale_fill_viridis_d(option = "turbo", drop = FALSE, guide = "none") +
  ggtitle("HCPYA") +
  theme_gray(base_size = 7) +
  ylab(NULL) +
  xlab("Proportion Participants Excluded") +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5), labels = c(0, 0.5))

b <- lost2 |>
  filter(dataset == "hcpd") |>
  ggplot(aes(y = task, x = lost)) +
  facet_wrap(~type) +
  geom_col(aes(fill = scan), position = position_dodge(preserve = "single")) +
  scale_fill_viridis_d(option = "turbo", drop = FALSE, guide = "none") +
  ggtitle("HCPD") +
  theme_gray(base_size = 7) +
  ylab(NULL) +
  xlab("Proportion Participants Excluded") +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5), labels = c(0, 0.5))

c <- lost2 |>
  filter(dataset == "hcpa") |>
  ggplot(aes(y = task, x = lost)) +
  facet_wrap(~type) +
  geom_col(aes(fill = scan), position = position_dodge(preserve = "single")) +
  scale_fill_viridis_d(option = "turbo", drop = FALSE, guide = "none") +
  ggtitle("HCPA") +
  theme_gray(base_size = 7) +
  ylab(NULL) +
  xlab("Proportion Participants Excluded") +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5), labels = c(0, 0.5))


d <- lost2 |>
  filter(dataset == "ukb") |>
  ggplot(aes(y = task, x = lost)) +
  facet_grid(ses ~ type) +
  geom_col(aes(fill = scan), position = position_dodge(preserve = "single")) +
  scale_fill_viridis_d(option = "turbo", drop = FALSE, guide = "none") +
  ggtitle("UKB") +
  theme_gray(base_size = 7) +
  ylab(NULL) +
  xlab("Proportion Participants Excluded") +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5), labels = c(0, 0.5))

e <- lost2 |>
  filter(dataset == "abcd") |>
  ggplot(aes(y = task, x = lost)) +
  facet_grid(ses ~ type) +
  geom_col(aes(fill = scan), position = position_dodge(preserve = "single")) +
  scale_fill_viridis_d(option = "turbo", drop = FALSE, guide = "none") +
  ggtitle("ABCD") +
  theme_gray(base_size = 7) +
  ylab(NULL) +
  xlab("Proportion Participants Excluded") +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5), labels = c(0, 0.5))

f <- lost2 |>
  filter(dataset == "spacetop") |>
  ggplot(aes(y = task, x = lost)) +
  facet_grid(ses ~ type, scales = "free_y") +
  geom_col(aes(fill = scan), position = position_dodge(preserve = "single")) +
  scale_fill_viridis_d(option = "turbo", drop = FALSE) +
  ggtitle("Spacetop") +
  theme_gray(base_size = 7) +
  ylab(NULL) +
  xlab("Proportion Participants Excluded") +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5), labels = c(0, 0.5)) +
  guides(fill = guide_legend(nrow = 1))


p <- a + b + c + d + e + f +
  plot_layout(
    guides = "collect",
    design = "
    555111
    555111
    222333
    666444
    666444
    "
  ) &
  theme(legend.position = "bottom")

ggsave(
  "figures/all_motion_exclusion2_filtered.png",
  p,
  device = ragg::agg_png,
  width = 6.5,
  height = 6.5
)
