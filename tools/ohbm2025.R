library(dplyr)
library(ggplot2)
library(tidyr)
library(ggdist)
library(patchwork)
library(purrr)

# plot_hcp_fd <- function(hcpya, hcp_design, task, fd_max_y=0.4){
#   hcpya |>
#     filter(stringr::str_detect(task, .env$task)) |>
#     ggplot(aes(x=t, y=median, color=ped)) +
#     facet_grid(ped~task, scales = "free_x") +
#     geom_rect(
#       data=filter(hcp_design, stringr::str_detect(task, .env$task)),
#       aes(
#         xmin=onset,
#         xmax = onset+duration,
#         ymin = Inf,
#         ymax = -Inf,
#         fill=type),
#       alpha = 0.5,
#       inherit.aes = FALSE) +
#     geom_line() +
#     ylab("Framewise Displacement") +
#     xlab("Volume") +
#     scale_fill_viridis_d(option = "turbo") +
#     theme(legend.position = "bottom") +
#     guides(color=FALSE) +
#     coord_cartesian(ylim = c(0, fd_max_y))
# }
#
# targets::tar_load(c(hcp_design, ukb_design, ukb_responses))
#
# hcpya <- arrow::open_dataset("data/motion/human-connectome-project-openaccess.parquet") |>
#   collect() |>
#   summarise(
#     median = mean(framewise_displacement),
#     .by = c(t, task, ped, ses)
#   ) |>
#   mutate(ses = glue::glue("Ses {ses}"))
#
#
# a <- plot_hcp_fd(hcpya, hcp_design, "WM") + ggtitle("Ses 1")
# b <- plot_hcp_fd(hcpya, hcp_design, "GAMBLING")
# c <- plot_hcp_fd(hcpya, hcp_design, "MOTOR")
#
# d <- plot_hcp_fd(hcpya, hcp_design, "LANGUAGE") + ggtitle("Ses 2")
# e <- plot_hcp_fd(hcpya, hcp_design, "SOCIAL")
# f <- plot_hcp_fd(hcpya, hcp_design, "RELATIONAL")
# g <- plot_hcp_fd(hcpya, hcp_design, "EMOTION")
#
# hcpya_p <- a + b + c + d + e + f + g + plot_layout(
#   design = "
#   123#
#   4567
#   ") &
#   theme(legend.position = "bottom")
#
# hcpya_p <- hcpya |>
#   filter(t>0) |>
#   ggplot() +
#   facet_wrap(~task, scales = "free_x") +
#   geom_line(aes(x=t, y=median, color=ped)) +
#   ylab("Framewise Displacement") +
#   xlab("Volume") +
#   scale_fill_viridis_d(option = "turbo") +
#   theme(legend.position = "bottom") +
#   coord_cartesian(ylim = c(0, 0.4))
#
#
# ukb <- arrow::open_dataset("data/motion/ukb.parquet") |>
#   filter(t>0) |>
#   filter(datatype=="20227" | t<=332) |>
#   mutate(task = if_else(datatype=="20227" ,"rest", "faces/shapes")) |>
#   collect() |>
#   summarise(
#     sem = mad(framewise_displacement) / sqrt(n()),
#     median = median(framewise_displacement),
#     .by = c(t, task, ses)
#   ) |>
#   mutate(ses = glue::glue("Ses {ses}"))
#
#
# ukb_p <- ukb |>
#   ggplot() +
#   facet_wrap(~task, scales = "free_x", ncol = 1) +
#   geom_line(aes(x=t, y=median, color=ses)) +
#   ylab("Framewise Displacement") +
#   xlab("Volume") +
#   scale_fill_viridis_d(option = "turbo") +
#   theme(legend.position = "bottom") +
#   coord_cartesian(ylim = c(0, 0.4))
#
#
# abcd_n <- arrow::open_dataset("data/motion/abcd.parquet") |>
#   summarise(n_tr = max(t), .by=c(sub, task, ses, run)) |>
#   count(n_tr, task, ses, run) |>
#   collect() |>
#   slice_max(order_by=n, n=1, by=c(task, ses, run)) |>
#   distinct(n_tr, task, ses, run) |>
#   arrow::as_arrow_table() |>
#   mutate(
#     ses=arrow::cast(ses, arrow::large_utf8()),
#     task=arrow::cast(task, arrow::large_utf8()),
#     run=arrow::cast(run, arrow::large_utf8())
#   )
#
# abcd <- arrow::open_dataset("data/motion/abcd.parquet") |>
#   left_join(abcd_n) |>
#   filter(t<n_tr, t>0) |>
#   collect() |>
#   summarise(
#     sem = mad(framewise_displacement) / sqrt(n()),
#     median = median(framewise_displacement),
#     .by = c(t, task, ses, run)
#   )
#
# abcd_p <- abcd |>
#   mutate(
#     ses=factor(ses, levels=c("baselineYear1Arm1","2YearFollowUpYArm1","4YearFollowUpYArm1")),
#   ) |>
#   filter(run %in% c("01","02","03","04")) |>
#   ggplot() +
#   facet_grid(task~ses, scales = "free_x") +
#   geom_line(aes(x=t, y=median, color=run)) +
#   ylab("Framewise Displacement") +
#   xlab("Volume") +
#   scale_fill_viridis_d(option = "turbo") +
#   theme(legend.position = "bottom") +
#   coord_cartesian(ylim = c(0, 0.4))
#
#
# hcpdev <- arrow::open_dataset("data/motion/HCPDevelopmentRec.parquet") |>
#   filter(t>0) |>
#   collect() |>
#   summarise(
#     sem = mad(framewise_displacement) / sqrt(n()),
#     median = median(framewise_displacement),
#     .by = c(t, task, ses, ped)
#   )
#
# hcpdev_p <- hcpdev  |>
#   ggplot() +
#   facet_wrap(~task, scales = "free_x") +
#   geom_line(aes(x=t, y=median, color=ped)) +
#   ylab("Framewise Displacement") +
#   xlab("Volume") +
#   scale_fill_viridis_d(option = "turbo") +
#   theme(legend.position = "bottom") +
#   coord_cartesian(ylim = c(0, 0.4))
#
#
# hcpaging <- arrow::open_dataset("data/motion/HCPAgingRec.parquet") |>
#   filter(t>0) |>
#   collect() |>
#   summarise(
#     sem = mad(framewise_displacement) / sqrt(n()),
#     median = median(framewise_displacement),
#     .by = c(t, task, ses, ped)
#   )
#
#
# hcpaging_p <- hcpaging  |>
#   ggplot() +
#   facet_wrap(~task, scales = "free_x") +
#   geom_line(aes(x=t, y=median, color=ped)) +
#   ylab("Framewise Displacement") +
#   xlab("Volume") +
#   scale_fill_viridis_d(option = "turbo") +
#   theme(legend.position = "bottom") +
#   coord_cartesian(ylim = c(0, 0.4))
#
#
# p <- (hcpya_p + ggtitle("HCPYA")) +
#   (hcpdev_p + ggtitle("HCPYD")) +
#   hcpaging_p + ggtitle("HCPA") +
#   (ukb_p + ggtitle("UKB")) +
#   (abcd_p + ggtitle("ABCD")) +
#   plot_layout(
#     guides = "collect",
#     design = "
#     111111
#     222333
#     444555
#     ") +
#   plot_annotation(tag_levels = "a") &
#   theme_gray(base_size = 7) &
#   theme(
#     legend.position = "bottom"
#   )
#
# ggsave(
#   "figures/all_motion.png",
#   device = ragg::agg_png,
#   width = 6.5,
#   height = 7)

targets::tar_load(by_time)

a <- by_time |>
  filter(dataset == "hcpya") |>
  ggplot() +
  facet_wrap(~task, scales = "free_x") +
  geom_line(aes(x = time, y = median, color = ped.run)) +
  ylab("Framewise Displacement") +
  xlab("Time (sec)") +
  coord_cartesian(ylim = c(0, 0.3)) +
  scale_colour_viridis_d(option = "turbo", drop = FALSE, guide = "none") +
  ggtitle("HCPYA")

b <- by_time |>
  filter(dataset == "hcpd") |>
  ggplot() +
  facet_wrap(~task, scales = "free_x") +
  geom_line(aes(x = time, y = median, color = ped.run)) +
  ylab("Framewise Displacement") +
  xlab("Time (sec)") +
  coord_cartesian(ylim = c(0, 0.3)) +
  scale_colour_viridis_d(option = "turbo", drop = FALSE, guide = "none") +
  ggtitle("HCPD")

c <- by_time |>
  filter(dataset == "hcpa") |>
  ggplot() +
  facet_wrap(~task, scales = "free_x") +
  geom_line(aes(x = time, y = median, color = ped.run)) +
  ylab("Framewise Displacement") +
  xlab("Time (sec)") +
  coord_cartesian(ylim = c(0, 0.3)) +
  scale_colour_viridis_d(option = "turbo", drop = FALSE, guide = "none") +
  ggtitle("HCPA")

d <- by_time |>
  filter(dataset == "ukb", ses == "2") |>
  ggplot() +
  facet_wrap(~task, scales = "free_x", nrow = 2) +
  geom_line(aes(x = time, y = median, color = ped.run)) +
  ylab("Framewise Displacement") +
  xlab("Time (sec)") +
  coord_cartesian(ylim = c(0, 0.3)) +
  scale_colour_viridis_d(option = "turbo", drop = FALSE, guide = "none") +
  ggtitle("UKB")

e <- by_time |>
  filter(dataset == "abcd") |>
  mutate(
    ses = forcats::fct_recode(
      ses,
      baseline = "baselineYear1Arm1",
      Year2 = "2YearFollowUpYArm1",
      Year4 = "4YearFollowUpYArm1"
    )
  ) |>
  bind_rows(
    tibble(
      ped.run = factor(c("1.1", "1.2", "2.3", "2.4")),
      ses = rep(factor("baseline"), times = 4),
      task = rep("rest", times = 4)
    )
  ) |>
  ggplot() +
  facet_grid(ses ~ task, scales = "free_x") +
  geom_line(aes(x = time, y = median, color = ped.run)) +
  ylab("Framewise Displacement") +
  xlab("Time (sec)") +
  coord_cartesian(ylim = c(0, 0.3)) +
  scale_colour_viridis_d(option = "turbo") +
  ggtitle("ABCD") +
  guides(color = guide_legend(nrow = 1))

p <- a + b + c + d + e +
  plot_layout(
    guides = "collect",
    design = "
    111111
    222333
    444555
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


targets::tar_load(c(by_run, demographics))

demographics <- demographics |>
  mutate(
    dataset = case_match(
      dataset,
      "HCPAgingRec" ~ "hcpa",
      "HCPDevelopmentRec" ~ "hcpd",
      "human_connectome_project_openaccess" ~ "hcpya",
      "UKB" ~ "ukb",
      "ABCD" ~ "abcd"
    )
  )


all_by_run <- by_run |>
  left_join(distinct(demographics, dataset, sub, ses, age, sex_gender, bmi)) |>
  filter(!ses == "baselineYear1Arm1") |>
  group_nest(dataset) |>
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
  mutate(
    dataset = case_match(
      dataset,
      "abcd" ~ "ABCD",
      "hcpa" ~ "HCPA",
      "hcpd" ~ "HCPD",
      "hcpya" ~ "HCPYA",
      "ukb" ~ "UKB"
    )
  )


a <- all_by_run |>
  filter(!is.na(age), !is.na(sex_gender)) |>
  mutate(N = n(), .by = c(dataset)) |>
  crossing(thresh = c(0.1, 0.2, 0.3, 0.4, 0.5, Inf)) |>
  summarise(
    n = unique(sum(loc > thresh) / N),
    .by = c(dataset, thresh, age)
  ) |>
  mutate(thresh = factor(thresh)) |>
  ggplot(aes(x = thresh, y = n)) +
  geom_line(aes(color = age, group = age)) +
  facet_wrap(~dataset, nrow = 1) +
  scale_color_viridis_d(option = "turbo", name = "Age Quartile") +
  xlab("Avg. Framewise Displacement Threshold") +
  ylab("Proportion Lost to Thresholding") +
  theme_gray(base_size = 9)

b <- all_by_run |>
  filter(!is.na(bmi), !is.na(sex_gender)) |>
  mutate(
    bmi = cut(
      bmi,
      breaks = c(0, 18.5, 25, 30, Inf),
      labels = c("Underweight", "Normal", "Overweight", "Obese"),
      ordered_result = TRUE
    )
  ) |>
  mutate(N = n(), .by = c(dataset)) |>
  crossing(thresh = c(0.1, 0.2, 0.3, 0.4, 0.5, Inf)) |>
  summarise(
    n = unique(sum(loc > thresh) / N),
    .by = c(dataset, thresh, bmi)
  ) |>
  mutate(thresh = factor(thresh)) |>
  ggplot(aes(x = thresh, y = n)) +
  geom_line(aes(color = bmi, group = bmi)) +
  facet_wrap(~dataset, nrow = 1) +
  scale_color_viridis_d(option = "turbo", name = "BMI") +
  xlab("Avg. Framewise Displacement Threshold") +
  ylab("Proportion Lost to Thresholding") +
  theme_gray(base_size = 9)

p <- a + b +
  plot_layout(nrow = 2) +
  plot_annotation(tag_levels = "a", tag_suffix = ")") &
  theme(legend.position = "bottom")


ggsave("figures/ohbm2.png", p, device = ragg::agg_png, height = 6, width = 6)


lost <- targets::tar_read(lost) |>
  mutate(
    ses = forcats::fct_recode(
      ses,
      baseline = "baselineYear1Arm1",
      Year2 = "2YearFollowUpYArm1",
      Year4 = "4YearFollowUpYArm1",
      "Ses 2" = "2",
      "Ses 3" = "3"
    ) |>
      forcats::fct_relevel("Year2", after = Inf) |>
      forcats::fct_relevel("Year4", after = Inf),
    scan = factor(scan)
  )

a <- lost |>
  filter(dataset == "hcpya") |>
  ggplot(aes(y = task, x = lost)) +
  facet_wrap(~type) +
  geom_col(aes(fill = scan), position = position_dodge(preserve = "single")) +
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
  scale_fill_viridis_d(option = "turbo", drop = FALSE) +
  ggtitle("ABCD") +
  theme_gray(base_size = 7) +
  ylab(NULL) +
  xlab("Proportion Participants Excluded") +
  xlim(0, 1)

p <- a + b + c + d + e +
  plot_layout(guides = "collect", nrow = 3) &
  theme(legend.position = "bottom")

ggsave(
  "figures/all_motion_exclusion.png",
  p,
  device = ragg::agg_png,
  width = 6.5,
  height = 6.5
)

targets::tar_load(c(hcpya, hcpa, hcpd, abcd, ukb))

by_prop <- purrr::map(
  list(
    hcpya = dplyr::mutate(hcpya, ses = "1"),
    hcpa = hcpa,
    hcpd = hcpd,
    ukb = ukb,
    abcd = abcd
  ),
  .summarize_by_prop,
  threshold = 0.2
) |>
  dplyr::bind_rows(.id = "dataset") |>
  dplyr::filter(lost > prop_thresh) |>
  dplyr::distinct(dataset, task, ses, run, ped, sub)

by_max <- purrr::map(
  list(
    hcpya = dplyr::mutate(hcpya, ses = "1"),
    hcpa = hcpa,
    hcpd = hcpd,
    ukb = ukb,
    abcd = abcd
  ),
  .summarize_by_max
) |>
  dplyr::bind_rows(.id = "dataset") |>
  dplyr::filter(max_fd > max_fd_thresh) |>
  dplyr::distinct(dataset, task, ses, run, ped, sub)

by_avg <- by_run |>
  dplyr::filter(loc > mfd_thresh) |>
  dplyr::distinct(dataset, task, ses, run, ped, sub)

n_subs <- by_run |>
  dplyr::count(dataset, task, ses, run, ped, name = "n_sub")

tmp <- dplyr::bind_rows(by_max, by_avg, by_prop) |>
  dplyr::distinct() |>
  mutate(lost = TRUE) |>
  right_join(distinct(demographics, dataset, sub, ses, age, sex_gender, bmi)) |>
  mutate(lost = if_else(is.na(lost), FALSE, lost)) |>
  filter(!is.na(dataset))



for_cor <- by_run |>
  left_join(distinct(demographics, dataset, sub, ses, age, sex_gender, bmi)) |>
  group_nest(dataset) |>
  mutate(data = map(data, ~ cor.test(.x$age, .x$loc) |> broom::tidy())) |>
  unnest(data)
