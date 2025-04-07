library(dplyr)
library(ggplot2)
library(tidyr)
library(ggdist)
library(patchwork)
library(scattermore)


counts <- arrow::open_dataset("~/data/motion/derivatives/abcd.parquet") |>
  group_by(sub, task, run, ses) |>
  summarise(
    N = n(),
    .groups = "drop"
  ) |>
  collect() |>
  count(ses, task, run, N) |>
  mutate(mode = n == max(n), .by = c(ses, task, run)) |>
  filter(mode) |>
  select(ses, task, run, n_t_expected = N) |>
  arrow::arrow_table() |>
  mutate(
    task = arrow:::cast(task, arrow::large_utf8()),
    ses = arrow:::cast(ses, arrow::large_utf8()),
    run = arrow:::cast(run, arrow::large_utf8()),
  )

image03 <- readr::read_tsv(
  "/Users/psadil/data/events/sourcedata/image03.tsv",
  col_select = c(src_subject_id, starts_with("scanner"))
) |>
  rename(sub = src_subject_id) |>
  mutate(sub = stringr::str_remove(sub, "_")) |>
  filter(!scanner_manufacturer_pd == "GE MEDICAL SYSTEMS") |>
  select(sub) |>
  arrow::as_arrow_table() |>
  transmute(sub = cast(sub, arrow::large_utf8()))

abcd <- arrow::open_dataset("~/data/motion/derivatives/abcd.parquet") |>
  left_join(counts) |>
  filter(t <= n_t_expected) |>
  filter(!run == "05", !run == "06") |> # just a few resting state scans
  group_by(t, task, run, ses, n_t_expected) |>
  summarise(
    framewise_displacement = mean(framewise_displacement),
    N = n(),
    .groups = "drop"
  ) |>
  collect()

abcd_events <- get_abcd_events()


p <- abcd |>
  mutate(
    ses = factor(
      ses,
      levels = c("baselineYear1Arm1", "2YearFollowUpYArm1", "4YearFollowUpYArm1"),
      ordered = TRUE
    ),
    task = factor(task, levels = c("rest", "sst", "mid", "nback"))
  ) |>
  filter(t > 0) |>
  mutate(t = t * 0.8) |>
  ggplot(aes(x = t, y = framewise_displacement, color = run)) +
  facet_grid(task ~ ses) +
  geom_rect(
    data = mutate(abcd_events,
      ses = factor(
        ses,
        levels = c("baselineYear1Arm1", "2YearFollowUpYArm1", "4YearFollowUpYArm1"),
        ordered = TRUE
      ),
      task = factor(task, levels = c("rest", "sst", "mid", "nback"))
    ),
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
  geom_line(linewidth = 1) +
  scale_fill_viridis_d(option = "turbo") +
  ylim(0, 1) +
  theme(legend.position = "bottom")

ggsave("figures/abcd_motion_events.png", p, width = 16)
