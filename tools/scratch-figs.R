library(dplyr)
library(ggplot2)
library(tidyr)
library(ggdist)
library(patchwork)
library(purrr)
library(gtsummary)

GeomSplitViolin <- ggproto(
  "GeomSplitViolin",
  GeomViolin,
  draw_group = function(self, data, ..., draw_quantiles = NULL) {
    data <- transform(data,
      xminv = x - violinwidth * (x - xmin),
      xmaxv = x + violinwidth * (xmax - x)
    )
    grp <- data[1, "group"]
    newdata <- plyr::arrange(
      transform(data, x = if (grp %% 2 == 1) xminv else xmaxv),
      if (grp %% 2 == 1) y else -y
    )
    newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
    newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])
    if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
      stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <= 1))
      quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
      aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
      aesthetics$alpha <- rep(1, nrow(quantiles))
      both <- cbind(quantiles, aesthetics)
      quantile_grob <- GeomPath$draw_panel(both, ...)
      ggplot2:::ggname(
        "geom_split_violin",
        grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob)
      )
    } else {
      ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
    }
  }
)

geom_split_violin <- function(mapping = NULL,
                              data = NULL,
                              stat = "ydensity",
                              position = "identity", ...,
                              draw_quantiles = NULL,
                              trim = TRUE,
                              scale = "area",
                              na.rm = FALSE,
                              show.legend = NA,
                              inherit.aes = TRUE) {
  layer(
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
      na.rm = na.rm, ...
    )
  )
}

targets::tar_load(c(by_run, demographics))

mriqc <- arrow::open_dataset("data/bold") |>
  select(id = `_id`, fd_mean) |>
  collect() |>
  mutate(dataset = "mriqc")

all_by_run <- by_run |>
  filter(filtered) |>
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
  scale_color_viridis_d(option = "turbo") +
  ylab("Average Framewise\nDisplacement (mm)") +
  xlab("Age (years)") +
  theme_gray(base_size = 11) +
  guides(color = guide_legend(nrow = 1))

ranges_bmi <- all_by_run |>
  filter(stringr::str_detect(dataset, "SpaceTop", TRUE)) |>
  summarise(
    xmin = quantile(bmi, 0.975, na.rm = TRUE),
    xmax = quantile(bmi, 0.025, na.rm = TRUE),
    .by = dataset
  ) |>
  left_join(select(ranges, dataset, y), by = join_by(dataset))


b1 <- all_by_run |>
  filter(!is.na(bmi)) |>
  filter(stringr::str_detect(dataset, "ABCD", TRUE)) |>
  ggplot(aes(x = bmi, color = dataset)) +
  scattermore::geom_scattermore(aes(y = loc), alpha = 0.5, pointsize = 1) +
  geom_segment(aes(x = xmin, xend = xmax, y = y), data = filter(ranges_bmi, !dataset == "ABCD")) +
  geom_smooth(aes(y = loc, group = dataset), method = "lm") +
  coord_cartesian(ylim = c(-0.1, 1)) +
  scale_color_viridis_d(option = "turbo", guide = "none", drop = FALSE) +
  ylab("Average Framewise\nDisplacement (mm)") +
  xlab("Body Mass Index") +
  theme_gray(base_size = 11)

b2 <- all_by_run |>
  filter(!is.na(bmi)) |>
  filter(stringr::str_detect(dataset, "ABCD")) |>
  ggplot(aes(x = bmi, color = dataset)) +
  scattermore::geom_scattermore(aes(y = loc), alpha = 0.5, pointsize = 1) +
  geom_smooth(aes(y = loc, group = dataset), method = "lm") +
  geom_segment(aes(x = xmin, xend = xmax, y = y), data = filter(ranges_bmi, dataset == "ABCD")) +
  coord_cartesian(ylim = c(-0.1, 1)) +
  scale_color_viridis_d(option = "turbo", guide = "none", drop = FALSE) +
  ylab("Average Framewise\nDisplacement (mm)") +
  xlab("Body Mass Index") +
  theme_gray(base_size = 11)

summary_data <- all_by_run |>
  filter(sex_gender %in% c("F", "M")) |>
  summarise(
    mean = median(loc),
    min = quantile(loc, 0.25),
    max = quantile(loc, 0.75),
    .by = c(sex_gender, dataset)
  )

c <- all_by_run |>
  filter(sex_gender %in% c("F", "M")) |>
  ggplot(aes(x = dataset, color = sex_gender, y = loc)) +
  geom_split_violin(trim = FALSE, alpha = 0.5) +
  geom_pointrange(
    data = summary_data,
    aes(dataset, mean, ymin = min, ymax = max, color = sex_gender),
    shape = 20, # 95,
    position = position_dodge(width = 0.25)
  ) +
  coord_cartesian(ylim = c(-0.1, 1)) +
  ylab("Average Framewise\nDisplacement (mm)") +
  xlab(NULL) +
  theme_gray(base_size = 11) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

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
  ylab("Average Framewise Displacement (mm)") +
  xlab(NULL) +
  theme_gray(base_size = 11)

p <- a + b + c + d +
  plot_annotation(tag_levels = "a") +
  plot_layout(nrow = 2, guides = "collect") &
  theme(legend.position = "bottom")

ggsave(
  "figures/age-chart.tex",
  p,
  tikzDevice::tikz,
  width = 6.5,
  height = 4
)


p <- a + b2 + b1 + c + d +
  plot_layout(nrow = 1, guides = "collect", widths = c(3, 1, 3, 3, 1)) &
  theme(legend.position = "bottom")
ggsave(
  "figures/age-chart.png",
  p,
  width = 13,
  height = 6,
  device = ragg::agg_png
)


demo_tbl <- all_by_run |>
  mutate(
    dataset = interaction(dataset, ses, lex.order = TRUE, drop = TRUE)
  ) |>
  distinct(dataset, age, sex_gender, bmi, ses, sub)

demo_tbl |>
  na.omit() |>
  filter(stringr::str_detect(dataset, "ABCD")) |>
  mutate(ses = forcats::fct_drop(ses)) |>
  select(age, sex_gender, bmi, ses) |>
  tbl_summary(by = ses, missing = "no") |>
  as_gt() |>
  gt::gtsave("figures/demographics-abcd.tex")

demo_tbl |>
  na.omit() |>
  filter(stringr::str_detect(dataset, "HCPD")) |>
  select(age, sex_gender, bmi, ses) |>
  mutate(ses = forcats::fct_drop(ses)) |>
  tbl_summary(by = ses, missing = "no") |>
  as_gt() |>
  gt::gtsave("figures/demographics-hcpd.tex")

demo_tbl |>
  na.omit() |>
  filter(stringr::str_detect(dataset, "HCPA")) |>
  select(age, sex_gender, bmi, ses) |>
  mutate(ses = forcats::fct_drop(ses)) |>
  tbl_summary(by = ses, missing = "no") |>
  as_gt() |>
  gt::gtsave("figures/demographics-hcpa.tex")

demo_tbl |>
  na.omit() |>
  filter(stringr::str_detect(dataset, "HCPYA")) |>
  select(age, sex_gender, bmi, ses) |>
  mutate(ses = forcats::fct_drop(ses)) |>
  tbl_summary(by = ses, missing = "no") |>
  as_gt() |>
  gt::gtsave("figures/demographics-hcpya.tex")

demo_tbl |>
  na.omit() |>
  filter(stringr::str_detect(dataset, "UKB")) |>
  select(age, sex_gender, bmi, ses) |>
  mutate(ses = forcats::fct_drop(ses)) |>
  tbl_summary(by = ses, missing = "no") |>
  as_gt() |>
  gt::gtsave("figures/demographics-ukb.tex")


demo_tbl |>
  filter(stringr::str_detect(dataset, "SpaceTop")) |>
  select(age, sex_gender, ses) |>
  na.omit() |>
  mutate(ses = forcats::fct_drop(ses)) |>
  tbl_summary(by = ses, missing = "no") |>
  as_gt() |>
  gt::gtsave("figures/demographics-spacetop.tex")


a <- by_run |>
  filter(dataset == "hcpya", !filtered) |>
  ggplot(aes(y = task, color = scan, x = loc)) +
  geom_boxplot(outliers = FALSE, position = position_dodge(preserve = "single")) +
  xlab("Framewise Displacement (mm)") +
  ylab(NULL) +
  coord_cartesian(xlim = c(0, 1)) +
  scale_colour_viridis_d(option = "turbo", drop = FALSE, guide = "none") +
  ggtitle("HCPYA")

b <- by_run |>
  filter(!is.na(scan), !filtered) |>
  filter(dataset == "hcpd") |>
  ggplot(aes(y = task, color = scan, x = loc)) +
  geom_boxplot(outliers = FALSE, position = position_dodge(preserve = "single")) +
  xlab("Framewise Displacement (mm)") +
  ylab(NULL) +
  coord_cartesian(xlim = c(0, 1)) +
  scale_colour_viridis_d(option = "turbo", drop = FALSE, guide = "none") +
  ggtitle("HCPD")

c <- by_run |>
  filter(dataset == "hcpa", !filtered) |>
  ggplot(aes(y = task, color = scan, x = loc)) +
  geom_boxplot(outliers = FALSE, position = position_dodge(preserve = "single")) +
  xlab("Framewise Displacement (mm)") +
  ylab(NULL) +
  coord_cartesian(xlim = c(0, 1)) +
  scale_colour_viridis_d(option = "turbo", drop = FALSE, guide = "none") +
  ggtitle("HCPA")

d <- by_run |>
  filter(dataset == "ukb", !filtered) |>
  ggplot(aes(y = task, color = scan, x = loc)) +
  facet_wrap(~ses, nrow = 2, labeller = "label_both") +
  geom_boxplot(outliers = FALSE, position = position_dodge(preserve = "single")) +
  xlab("Framewise Displacement (mm)") +
  ylab(NULL) +
  coord_cartesian(xlim = c(0, 1)) +
  scale_colour_viridis_d(option = "turbo", drop = FALSE, guide = "none") +
  ggtitle("UKB")

e <- by_run |>
  filter(!is.na(scan), !filtered) |>
  filter(dataset == "abcd") |>
  ggplot(aes(y = task, color = scan, x = loc)) +
  facet_wrap(~ses, nrow = 3, labeller = "label_both") +
  geom_boxplot(outliers = FALSE, position = position_dodge(preserve = "single")) +
  xlab("Framewise Displacement (mm)") +
  ylab(NULL) +
  coord_cartesian(xlim = c(0, 1)) +
  scale_colour_viridis_d(option = "turbo", drop = FALSE, guide = "none") +
  ggtitle("ABCD")

f <- by_run |>
  filter(!is.na(scan), !filtered) |>
  filter(dataset == "spacetop") |>
  ggplot(aes(y = task, color = scan, x = loc)) +
  facet_wrap(~ses, scales = "free_y", labeller = "label_both") +
  geom_boxplot(outliers = FALSE, position = position_dodge(preserve = "single")) +
  xlab("Framewise Displacement (mm)") +
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
    555544
    223311
    666666
    666666
    "
  ) +
  plot_annotation() &
  theme_gray(base_size = 10) &
  theme(
    legend.position = "bottom"
  )


ggsave(
  "figures/all_motion-avg.tex",
  p,
  tikzDevice::tikz,
  width = 6.5,
  height = 7.7
)


p <- by_run |>
  filter(!is.na(scan)) |>
  ggplot(aes(y = task, color = scan, x = loc)) +
  ggh4x::facet_nested(dataset + ses ~ filtered, scales = "free_y", labeller = "label_both") +
  geom_boxplot(outliers = FALSE, position = position_dodge(preserve = "single")) +
  xlab("Framewise Displacement (mm)") +
  ylab(NULL) +
  coord_cartesian(xlim = c(0, 1)) +
  scale_colour_viridis_d(option = "turbo") +
  guides(color = guide_legend(nrow = 1)) +
  theme_gray(base_size = 10) +
  theme(
    legend.position = "bottom"
  )

ggsave(
  "figures/all_motion-avg.png",
  p,
  device = ragg::agg_png,
  width = 6.5,
  height = 13
)


targets::tar_load(c(by_time, abcd_events, hcpya_events))

hcpya_events <- hcpya_events |>
  mutate(
    task = forcats::fct_relabel(task, stringr::str_to_lower),
    ymin = if_else(ped == "LR", 0, 0.15),
    ymax = if_else(ped == "LR", 0.15, 0.3)
  ) |>
  filter(stringr::str_detect(task, "language", TRUE))

a <- by_time |>
  filter(dataset == "hcpya", !filtered) |>
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
    inherit.aes = FALSE
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
  "figures/all_motion.png",
  p,
  device = ragg::agg_png,
  width = 6.5,
  height = 7.5
)

ggsave(
  "figures/all_motion.tex",
  p,
  tikzDevice::tikz,
  width = 6.5,
  height = 7.5
)

p <- by_time |>
  filter(dataset == "abcd") |>
  ggplot(aes(x = time, y = median)) +
  facet_grid(ses ~ task, scales = "free_x") +
  scale_y_continuous(
    "Framewise Displacement (mm)",
    limits = c(0, 0.3)
  ) +
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
  geom_line(aes(color = scan, linetype = filtered)) +
  xlab("Time (seconds)") +
  ggtitle("ABCD") +
  theme(
    legend.position = "bottom"
  )

ggsave(
  "figures/abcd_bytime.png",
  p,
  device = ragg::agg_png,
  width = 13,
  height = 6
)
targets::tar_load(ukb_events)
ukb_responses <- get_ukb_responses() |>
  mutate(task = "faces/shapes")
p <- by_time |>
  filter(dataset == "ukb") |>
  ggplot(aes(x = time, y = median)) +
  facet_grid(ses ~ task, scales = "free") +
  scale_y_continuous(
    "Framewise Displacement (mm)",
    limits = c(0, 0.3)
  ) +
  geom_rect(
    aes(
      xmin = onset,
      xmax = onset + duration,
      ymin = Inf,
      ymax = -Inf,
      fill = type
    ),
    data = mutate(ukb_events, task = "faces/shapes"),
    alpha = 0.5,
    inherit.aes = FALSE
  ) +
  geom_vline(data = ukb_responses, aes(xintercept = onset)) +
  geom_line(aes(color = scan, linetype = filtered)) +
  xlab("Time (seconds)") +
  ggtitle("UKB") +
  theme(
    legend.position = "bottom"
  )

ggsave(
  "figures/ukb_bytime.png",
  p,
  device = ragg::agg_png,
  width = 13,
  height = 6
)

targets::tar_load(lost)

p <- lost |>
  ggplot(aes(y = task, fill = scan, x = lost)) +
  ggh4x::facet_nested(dataset + ses ~ filtered + type, scales = "free_y", labeller = "label_both") +
  geom_col(aes(fill = scan), position = position_dodge(preserve = "single")) +
  geom_errorbarh(
    aes(xmin = lower, xmax = upper),
    position = position_dodge(preserve = "single"),
  ) +
  scale_fill_viridis_d(option = "turbo") +
  scale_x_continuous(NULL, breaks = c(0, 0.5), labels = c(0, 0.5)) +
  guides(fill = guide_legend(nrow = 1)) +
  theme_gray(base_size = 10) +
  theme(
    legend.position = "bottom"
  )

ggsave(
  "figures/all_motion_exclusion.png",
  p,
  device = ragg::agg_png,
  width = 6.5,
  height = 13
)


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
  "figures/all_motion_exclusion.png",
  p,
  device = ragg::agg_png,
  width = 6.5,
  height = 6.5
)

ggsave(
  "figures/all_motion_exclusion.tex",
  p,
  tikzDevice::tikz,
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

p <- lost2 |>
  ggplot(aes(y = task, fill = scan, x = lost)) +
  ggh4x::facet_nested(dataset + ses ~ type, scales = "free_y", labeller = "label_both") +
  geom_col(aes(fill = scan), position = position_dodge(preserve = "single")) +
  scale_fill_viridis_d(option = "turbo") +
  ylab(NULL) +
  scale_x_continuous(breaks = c(0, 0.5), labels = c(0, 0.5)) +
  guides(fill = guide_legend(nrow = 1)) +
  theme_gray(base_size = 10) +
  theme(
    legend.position = "bottom"
  )

ggsave(
  "figures/all_motion_exclusion2.png",
  p,
  device = ragg::agg_png,
  width = 6.5,
  height = 13
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
  "figures/all_motion_exclusion2.png",
  p,
  device = ragg::agg_png,
  width = 6.5,
  height = 6.5
)


# spectrum

p <- targets::tar_read(hcpa_spectrum) |>
  mutate(
    task = stringr::str_to_lower(task)
  ) |>
  ggplot(aes(y = avg, x = freq)) +
  geom_raster(aes(fill = pxx)) +
  scale_fill_viridis_c(
    "Relative Power",
    option = "turbo",
    breaks = c(0, 0.5, 1),
    labels = c(0, 0.5, 1)
  ) +
  facet_grid(task + ped ~ param) +
  ylab("(<- lower avg FD) Participant (higher avg FD ->)") +
  xlab("Frequency (Hz)") +
  theme_gray(base_size = 12) +
  theme(
    legend.position = "bottom",
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank()
  )

ggsave(
  "figures/hcpa_spectrum.png",
  device = ragg::agg_png,
  width = 10,
  height = 6.5
)


p <- targets::tar_read(hcpd_spectrum) |>
  filter(ped == "PA") |>
  ggplot(aes(y = avg, x = freq)) +
  geom_raster(aes(fill = pxx)) +
  scale_fill_viridis_c(
    "Relative Power",
    option = "turbo",
    breaks = c(0, 0.5, 1),
    labels = c(0, 0.5, 1)
  ) +
  facet_grid(task + ped ~ param) +
  ylab("(<- lower avg FD) Participant (higher avg FD ->)") +
  xlab("Frequency (Hz)") +
  theme_gray(base_size = 12) +
  theme(
    legend.position = "bottom",
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank()
  )

ggsave(
  "figures/hcpd_spectrum.png",
  device = ragg::agg_png,
  width = 10,
  height = 6.5
)

p <- targets::tar_read(abcd_spectrum) |>
  filter(stringr::str_detect(param, "trans_y")) |>
  filter(run == "01") |>
  ggplot(aes(y = avg, x = freq)) +
  geom_raster(aes(fill = pxx)) +
  scale_fill_viridis_c(
    "Relative Power",
    option = "turbo",
    breaks = c(0, 0.5, 1),
    labels = c(0, 0.5, 1)
  ) +
  facet_grid(ses ~ task, scales = "free_y") +
  ylab("(<- lower avg FD) Participant (higher avg FD ->)") +
  xlab("Frequency (Hz)") +
  theme_gray(base_size = 12) +
  theme(
    legend.position = "bottom",
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank()
  )

ggsave(
  "figures/abcd_spectrum.png",
  device = ragg::agg_png,
  width = 10,
  height = 6.5
)


p <- targets::tar_read(hcpya_spectrum) |>
  filter(ped == "LR") |>
  mutate(task = stringr::str_to_lower(task)) |>
  ggplot(aes(y = avg, x = freq)) +
  geom_raster(aes(fill = pxx)) +
  scale_fill_viridis_c(
    "Relative Power",
    option = "turbo",
    breaks = c(0, 0.5, 1),
    labels = c(0, 0.5, 1)
  ) +
  facet_grid(task ~ param, scales = "free_y") +
  ylab("(<- lower avg FD) Participant (higher avg FD ->)") +
  xlab("Frequency (Hz)") +
  theme_gray(base_size = 10) +
  theme(
    legend.position = "bottom",
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank()
  )

ggsave(
  "figures/hcpya_spectrum.png",
  p,
  device = ragg::agg_png,
  width = 10,
  height = 6.5
)

p <- targets::tar_read(spacetop_spectrum) |>
  filter(stringr::str_detect(task, "align|narrati", TRUE)) |>
  mutate(task = stringr::str_to_lower(task)) |>
  ggplot(aes(y = avg, x = freq)) +
  geom_raster(aes(fill = pxx)) +
  scale_fill_viridis_c(
    "Relative Power",
    option = "turbo",
    breaks = c(0, 0.5, 1),
    labels = c(0, 0.5, 1)
  ) +
  ggh4x::facet_nested(ses + task ~ param, scales = "free_y") +
  ylab("(<- lower avg FD) Participant (higher avg FD ->)") +
  xlab("Frequency (Hz)") +
  theme_gray(base_size = 10) +
  theme(
    legend.position = "bottom",
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank()
  )

ggsave(
  "figures/spacetop_spectrum.png",
  p,
  device = ragg::agg_png,
  width = 10,
  height = 6.5
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
            breaks = c(-Inf, quantile(age, probs = c(0.25, 0.5, 0.75, 1), na.rm = TRUE)),
            labels = c("[0, 25)", "[25, 50)", "[50, 75)", "[75, 100]"),
            ordered_result = TRUE
          )
        )
    )
  ) |>
  unnest(data) |>
  filter(!is.na(age), !is.na(sex_gender)) |>
  mutate(N = n(), .by = c(dataset, ses, sex_gender)) |>
  crossing(thresh = c(seq(0.1, 1, by = 0.1))) |>
  summarise(
    n = unique(sum(loc > thresh) / N),
    .by = c(dataset, thresh, age, ses, sex_gender)
  ) |>
  ggplot(aes(x = thresh, y = n)) +
  geom_line(aes(color = age, linetype = sex_gender)) +
  facet_wrap(~ dataset + ses, nrow = 1) +
  scale_color_viridis_d(option = "turbo", name = "Age Quartile") +
  xlab("Avg. Framewise Displacement Threshold") +
  ylab("Proportion Lost to Thresholding") +
  scale_x_continuous(
    "Avg. Framewise Displacement Threshold",
    limits = c(0, 1),
    breaks = c(0, 0.5),
    labels = c(0, 0.5)
  ) +
  theme_gray(base_size = 14) +
  theme(
    legend.position = "bottom"
  )

ggsave(
  "figures/motion-figs-2.png",
  p,
  device = ragg::agg_png,
  width = 13,
  height = 6.5
)

p <- all_by_run |>
  mutate(ses = forcats::fct_collapse(ses, " " = c("1", "2", "3", "4"))) |>
  filter(!is.na(bmi), !is.na(sex_gender)) |>
  mutate(
    bmi = cut(
      bmi,
      breaks = c(0, 18.5, 25, 30, Inf),
      labels = c("Underweight", "Normal", "Overweight", "Obese"),
      ordered_result = TRUE
    )
  ) |>
  mutate(N = n(), .by = c(dataset, ses, sex_gender)) |>
  crossing(thresh = c(seq(0.1, 1, by = 0.1))) |>
  summarise(
    n = unique(sum(loc > thresh) / N),
    .by = c(dataset, thresh, bmi, ses, sex_gender)
  ) |>
  ggplot(aes(x = thresh, y = n)) +
  geom_line(aes(color = bmi, linetype = sex_gender)) +
  facet_wrap(~ dataset + ses, nrow = 1) +
  scale_color_viridis_d(option = "turbo", name = "BMI") +
  ylab("Proportion Lost to Thresholding") +
  scale_x_continuous(
    "Avg. Framewise Displacement Threshold",
    limits = c(0, 1),
    breaks = c(0, 0.5),
    labels = c(0, 0.5)
  ) +
  theme_gray(base_size = 14) +
  theme(
    legend.position = "bottom"
  )

ggsave(
  "figures/motion-figs-2b.png",
  p,
  device = ragg::agg_png,
  width = 13,
  height = 6.5
)


# notch params

hcpa <- arrow::open_dataset("/Users/psadil/data/motion/derivatives/hcpa_spectrum.parquet") |>
  filter(between(freq, 0.1, 0.6), param == "trans_y") |>
  collect() |>
  slice_max(order_by = pxx, n = 1, by = c(sub, ses, task, ped)) |>
  summarise(
    freq = median(freq),
    .by = c(sub, task, param)
  ) |>
  summarise(
    lower = median(freq),
    upper = quantile(freq, 0.75),
    .by = c(task)
  )
hcpd <- arrow::open_dataset("/Users/psadil/data/motion/derivatives/hcpd_spectrum.parquet") |>
  filter(between(freq, 0.1, 0.6), param == "trans_y") |>
  collect() |>
  slice_max(order_by = pxx, n = 1, by = c(sub, ses, task, ped)) |>
  summarise(
    freq = median(freq),
    .by = c(sub, task)
  ) |>
  summarise(
    lower = median(freq),
    upper = quantile(freq, 0.75),
    .by = c(task)
  )
hcpya <- arrow::open_dataset("/Users/psadil/data/motion/derivatives/hcpya_spectrum.parquet") |>
  filter(between(freq, 0.1, 0.6), param == "trans_x") |>
  collect() |>
  slice_max(order_by = pxx, n = 1, by = c(sub, ses, task, ped)) |>
  summarise(
    freq = median(freq),
    .by = c(sub, task)
  ) |>
  summarise(
    lower = median(freq),
    upper = quantile(freq, 0.75),
    .by = c(task)
  )

abcd <- arrow::open_dataset("/Users/psadil/data/motion/derivatives/abcd_spectrum.parquet") |>
  filter(between(freq, 0.1, 0.6), param == "trans_y") |>
  collect() |>
  slice_max(order_by = pxx, n = 1, by = c(sub, ses, task)) |>
  summarise(
    lower = median(freq),
    upper = quantile(freq, 0.75),
    .by = c(task, ses)
  )

ukb <- arrow::open_dataset("/Users/psadil/data/motion/derivatives/ukb_spectrum.parquet") |>
  filter(between(freq, 0.1, 0.6), param == "trans_y") |>
  collect() |>
  slice_max(order_by = pxx, n = 1, by = c(sub, ses, task)) |>
  summarise(
    freq = median(freq),
    .by = c(sub, task)
  ) |>
  summarise(
    lower = median(freq),
    upper = quantile(freq, 0.75),
    .by = c(task)
  )

spacetop <- arrow::read_parquet("/Users/psadil/data/motion/derivatives/spacetop_spectrum.parquet") |>
  filter(between(freq, 0.1, 0.6), param == "trans_y") |>
  slice_max(order_by = pxx, n = 1, by = c(sub, ses, task, run)) |>
  summarise(
    freq = median(freq),
    .by = c(sub, task)
  ) |>
  summarise(
    lower = median(freq),
    upper = quantile(freq, 0.75),
    .by = c(task)
  )

dplyr::bind_rows(
  list(
    abcd = abcd,
    hcpa = hcpa,
    hcpd = hcpd,
    hcpya = hcpya,
    ukb = ukb,
    spacetop = spacetop
  ),
  .id = "dataset"
) |>
  DT::datatable(filter = "top")

dplyr::bind_rows(
  list(
    abcd = abcd,
    hcpa = hcpa,
    hcpd = hcpd,
    hcpya = hcpya,
    ukb = ukb,
    spacetop = spacetop
  ),
  .id = "dataset"
) |>
  summarise(
    lower = median(lower),
    upper = median(upper),
    .by = c(dataset, ses)
  )
