library(dplyr)

subs <- duckplyr::read_parquet_duckdb("data/Schaefer7n100p.parquet") |>
  distinct(sub) |>
  collect()

ukb <- targets::tar_read(ukb) |>
  filter(ses == "2", task == "rest") |>
  select(sub, t, framewise_displacement) |>
  mutate(sub = as.numeric(sub)) |>
  semi_join(subs, by = join_by(sub))

timeseries <- duckplyr::read_parquet_duckdb("data/Schaefer7n100p.parquet") |>
  filter(sub == 1007637) |>
  collect() |>
  left_join(ukb)


out2 <- get_cor_by_thresh(filter(out, sub == 1001445) |> mutate(ses = "2"))

library(ggplot2)

targets::tar_load(qc_ukb_summary)

tmp <- qc_ukb_summary |>
  tidyr::pivot_longer(avg:v) |>
  mutate(filtered = if_else(filtered, "filtered", "unfiltered"))

orig <- tmp |>
  dplyr::filter(iter == 0)

tmp |>
  filter(iter > 0) |>
  ggplot(aes(x = max_fd, y = value)) +
  scattermore::geom_scattermore(pointsize = 2) +
  geom_line(data = orig, color = "gold") +
  facet_grid(name ~ type + filtered, scales = "free_y") +
  theme_gray(base_size = 9)

ggsave("figures/surrogate.png", height = 7, width = 7, device = ragg::agg_png)

library(ggplot2)
library(dplyr)
library(tidyr)

targets::tar_load(qc_ukb_summary)

tmp <- qc_ukb_summary |>
  tidyr::pivot_longer(avg:v) |>
  mutate(filtered = if_else(filtered, "filtered", "unfiltered")) |>
  filter(sub == 1000974)

orig <- tmp |>
  dplyr::filter(iter == 0)

tmp |>
  filter(iter > 0) |>
  ggplot(aes(x = max_fd, y = value)) +
  scattermore::geom_scattermore(pointsize = 5) +
  geom_line(data = orig, color = "gold") +
  facet_grid(name ~ type + filtered, scales = "free_y") +
  theme_gray(base_size = 9)


tmp <- qc_ukb_summary |>
  mutate(filtered = if_else(filtered, "filtered", "unfiltered"))

orig <- tmp |>
  dplyr::filter(iter == 0) |>
  mutate(
    max_fd = cut(max_fd, c(-0.1, seq(0, 1, length.out = 11), Inf)),
    .by = c(sub, type, filtered)
  )


tmp |>
  filter(iter > 0) |>
  mutate(
    max_fd = cut(max_fd, c(-0.1, seq(0, 1, length.out = 11), Inf)),
    .by = c(sub, type, filtered)
  ) |>
  group_nest(sub, max_fd, type, filtered) |>
  mutate(edf = purrr::map(data, ~ ecdf(.x$avg))) |>
  select(-data) |>
  left_join(orig) |>
  mutate(
    rank = purrr::map2_dbl(edf, avg, ~ .x(.y)) |>
      cut(c(-0.01, seq(0, 1, length.out = 21))),
  ) |>
  count(rank, max_fd, type, filtered) |>
  mutate(
    value = n / sum(n),
    .by = c(max_fd, type, filtered)
  ) |>
  select(-n) |>
  complete(
    rank,
    max_fd,
    type,
    filtered,
    fill = list(value = 0)
  ) |>
  ggplot(aes(x = max_fd, y = rank)) +
  geom_raster(aes(fill = value)) +
  facet_wrap(~ type + filtered) +
  scale_fill_viridis_c(option = "turbo")

# experiment with MAC
library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)

targets::tar_load(mac)

mac |>
  mutate(
    type = replace_values(
      type,
      "clean" ~ "fix+smooth",
      "raw" ~ "24M+4S+BP+smooth",
      "rawraw" ~ "BP+smooth"
    ),
    filtered = replace_values(
      filtered,
      "TRUE" ~ "FD-filtered",
      "FALSE" ~ "FD",
      "clean-DPD" ~ "DPD-cleaned",
      "raw-DPD" ~ "DPD",
    )
  ) |>
  filter(
    !type == "BP+smooth",
    !filtered == "DPD-cleaned",
    !filtered == "FD-filtered"
  ) |>
  ggplot(aes(x = threshold, y = mac, linetype = filtered)) +
  geom_ribbon(
    aes(ymin = mac - 2 * sem, ymax = mac + 2 * sem, fill = type),
    alpha = 0.25
  ) +
  geom_line(aes(color = type)) +
  xlab("Proportion Frames Removed") +
  ylab("Mean Absolute Change (Scrubbed-Random)")

ggsave("figures/mac.png", device = ragg::agg_png, width = 6, height = 3)

# targets::tar_load(c(mac0, by_run))

# mac0 |>
#   filter(!filtered) |>
#   pivot_wider(names_from = type, values_from = mac) |>
#   mutate(clean_raw = clean - raw) |>
#   left_join(
#     by_run |>
#       filter(dataset == "ukb", ses == "2", task == "rest") |>
#       mutate(sub = as.integer(sub))
#   ) |>
#   ggplot(aes(x = loc, y = clean_raw)) +
#   facet_wrap(~threshold) +
#   geom_point(alpha = 0.5) +
#   geom_smooth(method = "lm")
