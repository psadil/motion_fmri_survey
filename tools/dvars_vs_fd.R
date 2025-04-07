library(tidyr)
library(dplyr)
library(ggplot2)
library(purrr)
library(stringr)

get_overlap <- function(motion, dvars) {
  motion_by <- motion |>
    semi_join(select(dvars, sub, datatype)) |>
    mutate(
      N = n(),
      q = rank(framewise_displacement) / N,
      .by = datatype
    ) |>
    crossing(thresh = 1 - seq(0, 0.5, by = 0.05)) |>
    filter(q >= thresh) |>
    mutate(metric = "framewise_displacement") |>
    select(-framewise_displacement)

  dvars_by <- dvars |>
    semi_join(select(motion, sub, datatype)) |>
    mutate(
      N = n(),
      q = rank(DPD) / N,
      .by = datatype
    ) |>
    crossing(thresh = 1 - seq(0, 0.5, by = 0.05)) |>
    filter(q >= thresh) |>
    mutate(metric = "DPD") |>
    select(-DPD)

  bind_rows(dvars_by, motion_by) |>
    select(-N, -q) |>
    pivot_wider(names_from = metric, values_from = sub, values_fn = list) |>
    rowwise() |>
    mutate(
      one_only = length(symdiff(framewise_displacement, DPD)),
      both = length(intersect(framewise_displacement, DPD)),
    ) |>
    select(-DPD, -framewise_displacement) |>
    pivot_longer(c(one_only, both))
}

ukb_dvars <- arrow::open_dataset("data/dvars/dataset=ukb", format = "arrow") |>
  mutate(
    keep = (stringr::str_detect(src, "20249") & t < 332) |
      (stringr::str_detect(src, "20227") & t < 490)
  ) |>
  filter(ses == 2, t > 0, keep, !src == "20227-filtered_func_data") |>
  summarise(
    DPD = sqrt(sum(DPD^2)),
    .by = c(sub, src)
  ) |>
  collect() |>
  mutate(datatype = str_extract(src, "[[:digit:]]{5}")) |>
  select(-src)


ukb_motion <- arrow::open_dataset("/Users/psadil/data/motion/derivatives/ukb.parquet") |>
  mutate(
    keep = (stringr::str_detect(datatype, "20249") & t < 332) |
      (stringr::str_detect(datatype, "20227") & t < 490)
  ) |>
  filter(ses == 2, t > 0, keep) |>
  summarise(
    framewise_displacement = mean(framewise_displacement),
    .by = c(datatype, sub)
  ) |>
  collect() |>
  mutate(sub = as.integer(sub))



overlap <- get_overlap(dvars = ukb_dvars, motion = ukb_motion)

overlap |>
  mutate(value = if_else(datatype == 20249, value / 38348, value / 44311)) |>
  ggplot(aes(x = thresh, y = value, color = name)) +
  facet_wrap(~datatype) +
  geom_line(alpha = 0.5) +
  ylab("Proportion Excluded")




p <- overlap2 |>
  ggplot(aes(x = thresh, y = value, color = measure)) +
  facet_grid(n_thresh ~ datatype) +
  geom_line(alpha = 0.5) +
  ylab("Proportion Participants Excluded") +
  xlab("Framewise Quantile Threshold") +
  theme(
    legend.position = "bottom"
  )

ggsave("figures/ukb_fd-v-dvars.png", device = ragg::agg_png, height = 6)
