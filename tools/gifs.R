library(dplyr)
library(ggplot2)

d <- nio::to_tbl("data/sub-2820279_ses-3_task-rest_bold.nii.gz")

p <- d |>
  filter(k == 38) |>
  ggplot(aes(x = i, y = j, fill = value)) +
  geom_raster() +
  scale_fill_viridis_c(option = "turbo", guide = "none") +
  coord_fixed() +
  theme_void() +
  gganimate::transition_time(t) +
  ggtitle(
    "Frame {frame} of {nframes}"
  )

gganimate::animate(p, fps = 1.25 * 2, device = "ragg_png", nframes = 490, renderer = gganimate::gifski_renderer("star.gif"))

pi <- d |>
  filter(i == 31) |>
  ggplot(aes(x = j, y = k, fill = value)) +
  geom_raster() +
  scale_fill_viridis_c(option = "turbo", guide = "none") +
  coord_fixed() +
  theme_void() +
  gganimate::transition_time(t) +
  ggtitle(
    "Frame {frame} of {nframes}"
  )

gganimate::animate(pi, fps = 1.25 * 2, device = "ragg_png", nframes = 490, renderer = gganimate::gifski_renderer("star-i.gif"))

pj <- d |>
  filter(j == 31) |>
  ggplot(aes(x = i, y = k, fill = value)) +
  geom_raster() +
  scale_fill_viridis_c(option = "turbo", guide = "none") +
  coord_fixed() +
  theme_void() +
  gganimate::transition_time(t) +
  ggtitle(
    "Frame {frame} of {nframes}"
  )

gganimate::animate(pj, fps = 1.25 * 2, device = "ragg_png", nframes = 490, renderer = gganimate::gifski_renderer("star-j.gif"))


ukb <- targets::tar_read(ukb) |>
  filter(ses == "3") |>
  filter(sub == "2820279") |>
  filter(task == "rest")

p2 <- ukb |>
  ggplot(aes(x = time, y = framewise_displacement, group = "run")) +
  geom_line() +
  ylab("Framewise Displacement\n(mm)") +
  xlab("Time (sec)") +
  theme_gray(base_size = 14)

ggsave(
  "figures/star2.png",
  p2,
  device = ragg::agg_png,
  width = 12,
  height = 3
)

ukb2 <- targets::tar_read(ukb2) |>
  filter(ses == "3") |>
  filter(sub == "2820279") |>
  filter(task == "rest") |>
  select(time, framewise_displacement_filtered2 = framewise_displacement_filtered)


p3 <- ukb |>
  select(time, starts_with("framewise")) |>
  left_join(ukb2) |>
  tidyr::pivot_longer(starts_with("frame"), names_to = "filtered") |>
  mutate(
    filtered = case_match(
      filtered,
      "framewise_displacement" ~ "raw",
      "framewise_displacement_filtered" ~ "[50-75%]",
      "framewise_displacement_filtered2" ~ "[25-75%]"
    )
  ) |>
  ggplot(aes(x = time, y = value, color = filtered)) +
  geom_line(alpha = 0.5) +
  ylab("Framewise Displacement (mm)") +
  xlab("Time (seconds)") +
  theme_gray(base_size = 14) +
  theme(legend.position = "bottom") +
  geom_point(alpha = 0)

p3_marg <- ggExtra::ggMarginal(p3,
  groupColour = TRUE,
  margins = "y", groupFill = TRUE, type = "density"
)

ggsave(
  "figures/ukb_motion_ex.png",
  p3_marg,
  device = ragg::agg_png,
  width = 12,
  height = 6
)

d <- nio::to_tbl("data/sub-NDARINVFX83EZC7_ses-baselineYear1Arm1_task-rest_run-01_bold.nii.gz")

pi <- d |>
  filter(i == 45) |>
  ggplot(aes(x = j, y = k, fill = value)) +
  geom_raster() +
  scale_fill_viridis_c(option = "turbo", guide = "none") +
  coord_fixed() +
  theme_void() +
  gganimate::transition_time(t) +
  ggtitle(
    "Frame {frame} of {nframes}"
  )

gganimate::animate(
  pi,
  fps = (1 / .72) * 2,
  device = "ragg_png",
  nframes = n_distinct(d$t),
  renderer = gganimate::av_renderer("abcd-i.mp4")
)

abcd <- targets::tar_read(abcd) |>
  filter(sub == "NDARINVFX83EZC7") |>
  filter(task == "rest") |>
  filter(ses == "baselineYear1Arm1") |>
  filter(run == "1")

p2 <- abcd |>
  ggplot(aes(x = t, y = framewise_displacement)) +
  geom_line() +
  ylab("Framewise Displacement (mm)") +
  xlab("Frame") +
  theme_gray(base_size = 14)

ggsave("abcd2.png", p2, device = ragg::agg_png, width = 11, height = 3)
