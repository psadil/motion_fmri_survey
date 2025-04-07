library(dplyr)
library(tidyr)
library(ggplot2)

read_timeseries <- function(src) {
  arrow::read_ipc_file(src, col_select = -t) |>
    corrr::correlate(quiet = TRUE) |>
    corrr::shave() |>
    corrr::stretch() |>
    na.omit()
}

ukb <- targets::tar_read(by_run) |>
  filter(dataset == "ukb", filtered, task == "rest") |>
  select(-sem, -filtered, -scan, -dataset, -task)

d <- tibble::tibble(src = fs::dir_ls("~/Desktop/dcs07/derivatives", glob = "*arrow", recurse = TRUE)) |>
  mutate(
    data = purrr::map(src, read_timeseries)
  )



# la <- freesurferformats::read.fs.annot("~/Downloads/lh.Schaefer2018_100Parcels_7Networks_order.annot")
# l <- freesurferformats::read.fs.surface("~/Downloads/fsaverage.L.sphere.164k_fs_L.surf.gii")
# vertices <- l$vertices
# colnames(vertices) <- c("x", "y", "z")
#
# ra <- freesurferformats::read.fs.annot("~/Downloads/rh.Schaefer2018_100Parcels_7Networks_order.annot")
# r <- freesurferformats::read.fs.surface("~/Downloads/fsaverage.R.sphere.164k_fs_R.surf.gii")
# r_vertices <- r$vertices
# colnames(r_vertices) <- c("x", "y", "z")
#
# positions <- bind_rows(
#     mutate(tibble::as_tibble(vertices), label_name = la$label_names),
#     mutate(tibble::as_tibble(r_vertices), label_name = ra$label_names)
#   ) |>
#   summarise(across(c(x,y,z), mean), .by = label_name) |>
#   arrange(label_name)
# p <- as.matrix(positions[2:4])
# rownames(p) <- positions$label_name
# distances <- corrr::as_cordf(as.matrix(dist(p, upper=TRUE))) |>
#   corrr::shave() |>
#   corrr::stretch() |>
#   na.omit() |>
#   rename(d=r)

positions <- readr::read_csv("~/Downloads/Schaefer2018_100Parcels_7Networks_order_FSLMNI152_1mm.Centroid_RAS.csv")
p <- as.matrix(positions[3:5])
rownames(p) <- positions$`ROI Name`
distances <- corrr::as_cordf(as.matrix(dist(p))) |>
  corrr::shave() |>
  corrr::stretch() |>
  na.omit() |>
  rename(d = r)


out <- d |>
  mutate(
    sub = stringr::str_extract(src, "(?<=sub=)[[:digit:]]+"),
    ses = stringr::str_extract(src, "(?<=ses=)[[:digit:]]+")
  ) |>
  left_join(ukb) |>
  filter(!is.na(loc)) |>
  select(-src) |>
  tidyr::unnest(data) |>
  left_join(distances) |>
  left_join(distances, by = join_by(x == y, y == x)) |>
  mutate(d = if_else(is.na(d.x), d.y, d.x)) |>
  select(-d.x, -d.y)


out |>
  filter(loc > quantile(loc, 0.95)) |>
  summarise(q = cor(r, loc), .by = c(x, y, d)) |>
  ggplot(aes(x = d, y = q)) +
  geom_point()


tests <- out |>
  group_nest(x, y, d) |>
  mutate(fit = purrr::map(data, ~ cor.test(.x$r, .x$loc) |> broom::tidy())) |>
  select(-data) |>
  unnest(fit)

tests |>
  summarise(sig = mean(p.value < 0.05))

lm(estimate ~ d, data = tests)




ukb <- targets::tar_read(ukb) |>
  filter(ses == "2") |>
  filter(stringr::str_detect(task, "rest")) |>
  select(sub, ses, starts_with("frame"), t) |>
  pivot_longer(starts_with("frame"), names_to = "filtered", values_to = "framewise_displacement") |>
  mutate(filtered = stringr::str_detect(filtered, "filt", TRUE))

d <- tibble::tibble(src = fs::dir_ls("~/Desktop/dcs07/derivatives", glob = "*arrow", recurse = TRUE)) |>
  slice_head(n = 1) |>
  crossing(window_start = seq(0, 200 - 151)) |>
  mutate(
    data = purrr::map2(src, window_start, ~ get_cor_by_thresh(src = .x, fd = ukb, window_start = .y))
  )


baseline <- d$data[[1]]

d |>
  select(-src) |>
  unnest(data) |>
  left_join(baseline, by = join_by(x, y)) |>
  summarise(
    r = cor(r.x, r.y),
    .by = c(max_fd.x, window_start)
  ) |>
  mutate(r = r - max(r)) |>
  ggplot(aes(x = max_fd.x, y = r)) +
  geom_point()


qcs <- targets::tar_read(qcs)

orig <- qcs3 |>
  filter(sub %in% c("3347745", "3554250", "2622095", "3496630", "4585017")) |>
  filter(iter == 0)

p <- qcs3 |>
  filter(sub %in% c("3347745", "3554250", "2622095", "3496630", "4585017")) |>
  filter(iter > 0) |>
  ggplot(aes(x = max_fd, y = r)) +
  facet_grid(sub ~ filtered, labeller = label_both) +
  geom_point(alpha = 0.1) +
  geom_point(alpha = 1, data = orig, color = "goldenrod")

ggsave("figures/qcrsfc-five-subs.png", device = ragg::agg_png)

centers <- qcs3 |>
  filter(
    sub %in% c("3347745", "3554250", "2622095", "3496630", "4585017"),
    iter > 0
  ) |>
  summarise(
    center = mean(r),
    .by = c(window_start, filtered, sub, ses)
  )

qcs3 |>
  filter(sub %in% c("3347745", "3554250", "2622095", "3496630", "4585017")) |>
  left_join(centers) |>
  mutate(
    r = r - center,
    group = if_else(iter == 0, "raw", "surrogate"),
    alpha = if_else(iter > 0, 0.01, 1)
  ) |>
  ggplot(aes(x = max_fd, y = r, color = group)) +
  facet_grid(sub ~ filtered, labeller = label_both) +
  geom_point(aes(alpha = alpha))



mutate(
  group = if_else(iter == 0, "raw", "surrogate"),
  alpha = if_else(iter > 0, 0.01, 1)
) |>
  ggplot(aes(x = max_fd, y = r, color = group)) +
  facet_grid(sub ~ filtered, labeller = label_both) +
  geom_point(aes(alpha = alpha))


ref <- qcs |> filter(iter == 0)
sur <- qcs |>
  select(-window_start) |>
  filter(iter > 0) |>
  group_nest(sub, ses, max_fd, filtered) |>
  left_join(ref) |>
  mutate(
    q = purrr::map2_dbl(data, r, ~ ecdf(.x$r)(.y))
  ) |>
  select(-data) |>
  mutate(
    max_fd = cut(max_fd, breaks = c(-Inf, seq(0, 1, by = 0.05), Inf)),
    q = cut(q, breaks = seq(0, 1, by = 0.05))
  ) |>
  count(q, max_fd, filtered) |>
  mutate(p = n / sum(n), .by = c(max_fd, filtered)) |>
  na.omit()

max_fds <- unique(sur$max_fd)

sur |>
  ggplot(aes(x = max_fd, y = q, fill = p)) +
  facet_wrap(~filtered) +
  geom_raster() +
  scale_fill_viridis_c(option = "turbo") +
  scale_x_discrete(
    breaks = c("(0,0.05]", "(0.2,0.25]", "(0.4,0.45]", "(0.6,0.65]", "(0.8,0.85]"),
    labels = c("0", "0.2", "0.4", "0.6", "0.8")
  )

ggsave("figures/qcs-surrogate.png", device = ragg::agg_png)


targets::tar_load(qcs3)

ref3 <- qcs3 |> filter(iter == 0)
sur3 <- qcs3 |>
  select(-window_start) |>
  filter(iter > 0) |>
  group_nest(sub, ses, max_fd, filtered) |>
  left_join(ref3) |>
  mutate(
    q = purrr::map2_dbl(data, r, ~ ecdf(.x$r)(.y))
  ) |>
  select(-data) |>
  mutate(
    max_fd = cut(max_fd, breaks = c(-Inf, seq(0, 1, by = 0.05), Inf)),
    q = cut(q, breaks = seq(0, 1, by = 0.05))
  ) |>
  count(q, max_fd, filtered) |>
  mutate(p = n / sum(n), .by = c(max_fd, filtered)) |>
  na.omit()

sur3 |>
  ggplot(aes(x = max_fd, y = q, fill = p)) +
  facet_wrap(~filtered) +
  geom_raster() +
  scale_fill_viridis_c(option = "turbo") +
  scale_x_discrete(
    breaks = c("(0,0.05]", "(0.2,0.25]", "(0.4,0.45]", "(0.6,0.65]", "(0.8,0.85]"),
    labels = c("0", "0.2", "0.4", "0.6", "0.8")
  )

ggsave("figures/sur3.png", device = ragg::agg_png, width = 5, height = 4.5)


sur3 |>
  mutate(max_fd = as.integer(max_fd)) |>
  ggplot(aes(x = max_fd, y = p, color = filtered)) +
  facet_wrap(~q) +
  geom_line()




targets::tar_load(qcs)

ref <- qcs |> filter(iter == 0)
sur <- qcs |>
  select(-window_start) |>
  filter(iter > 0) |>
  group_nest(sub, ses, max_fd, filtered) |>
  left_join(ref) |>
  mutate(
    q = purrr::map2_dbl(data, r, ~ ecdf(.x$r)(.y))
  ) |>
  select(-data) |>
  mutate(
    max_fd = cut(max_fd, breaks = c(-Inf, seq(0, 1, by = 0.05), Inf)),
    q = cut(q, breaks = seq(0, 1, by = 0.05))
  ) |>
  count(q, max_fd, filtered) |>
  mutate(p = n / sum(n), .by = c(max_fd, filtered)) |>
  na.omit()

sur |>
  ggplot(aes(x = max_fd, y = q, fill = p)) +
  facet_wrap(~filtered) +
  geom_raster() +
  scale_fill_viridis_c(option = "turbo") +
  scale_x_discrete(
    breaks = c("(0,0.05]", "(0.2,0.25]", "(0.4,0.45]", "(0.6,0.65]", "(0.8,0.85]"),
    labels = c("0", "0.2", "0.4", "0.6", "0.8")
  )


ggsave("figures/sur-avg-diff.png", device = ragg::agg_png)

ukb <- targets::tar_read(ukb) |>
  select(sub, ses, task, t, framewise_displacement, framewise_displacement_filtered) |>
  mutate(
    ses = as.integer(ses),
    sub = as.integer(sub),
    src = case_match(
      task,
      "rest" ~ "20227-filtered_func_data_clean",
      "faces/shapes" ~ "20249-filtered_func_data"
    )
  )

dvars <- arrow::open_dataset("data/dvars/dataset=ukb/", format = "arrow") |>
  select(t, DPD, ZD, sub, src, ses) |>
  na.omit() |>
  collect() |>
  left_join(select(ukb, -task), by = join_by(sub, ses, src, t))

tmp <- dvars |>
  na.omit() |>
  summarise(
    fd = cor(framewise_displacement, ZD),
    fdf = cor(framewise_displacement_filtered, ZD),
    fd2 = cor(framewise_displacement, DPD),
    fdf2 = cor(framewise_displacement_filtered, DPD),
    .by = c(sub, src, ses)
  )

tmp |>
  filter(ses == 2) |>
  pivot_longer(c(fd, fdf, fd2, fdf2)) |>
  ggplot(aes(x = value, fill = name)) +
  geom_histogram(position = "dodge") +
  facet_wrap(~src)


tmp |>
  pivot_longer(c(fd, fdf, fd2, fdf2)) |>
  summarise(
    value = median(value),
    .by = c(src, name)
  )




# targets::tar_load(qcs3_hcpya)

qcs3_hcpya <- arrow::open_dataset(fs::dir_ls("data/qcrsfc_150_01/", glob = "*parquet")) |>
  dplyr::collect()

subs <- c("128329", "919966", "107220", "123117", "123420", "140925")
orig <- qcs3_hcpya |>
  filter(sub %in% subs) |>
  filter(iter == 0)

p <- qcs3_hcpya |>
  filter(sub %in% subs) |>
  filter(iter > 0) |>
  ggplot(aes(x = max_fd, y = r)) +
  facet_grid(filtered + cleaned ~ sub, labeller = label_both, scales = "free") +
  scattermore::geom_scattermore(pointsize = 2) +
  geom_line(data = orig, color = "goldenrod")

ggsave("figures/qcrsfc-five-subs.png", p, device = ragg::agg_png, width = 8, height = 6)

qcs3_hcpya |>
  filter(sub %in% subs) |>
  filter(iter > 0) |>
  ggplot(aes(x = mean_fd, y = r)) +
  facet_grid(filtered + cleaned ~ sub, labeller = label_both, scales = "free") +
  scattermore::geom_scattermore(pointsize = 2) +
  geom_line(data = orig, color = "goldenrod")


ref <- qcs3_hcpya |> filter(iter == 0)
sur <- qcs3_hcpya |>
  select(-window_start) |>
  filter(iter > 0) |>
  group_nest(sub, max_fd, filtered, cleaned) |>
  left_join(ref) |>
  mutate(
    q = purrr::map2_dbl(data, r, ~ ecdf(.x$r)(.y))
  ) |>
  select(-data) |>
  mutate(
    max_fd = cut(max_fd, breaks = c(seq(0, 1, by = 0.05), Inf)),
    q = cut(q, breaks = c(-0.05, seq(0.05, 0.95, by = 0.05), 1.05))
  ) |>
  count(q, max_fd, filtered, cleaned) |>
  mutate(p = n / sum(n), .by = c(max_fd, filtered, cleaned))

p <- sur |>
  ggplot(aes(x = max_fd, y = q, fill = p)) +
  facet_grid(cleaned ~ filtered, labeller = "label_both") +
  geom_raster() +
  scale_fill_viridis_c(option = "turbo") +
  scale_x_discrete(
    breaks = c("(0,0.05]", "(0.1,0.15]", "(0.2,0.25]", "(0.3,0.35]", "(0.4,0.45]", "(0.5,0.55]", "(0.6,0.65]", "(0.7,0.75]", "(0.8,0.85]", "(0.9,0.95]"),
    labels = c("0", "0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9")
  ) +
  scale_y_discrete(
    breaks = c("(-0.05,0.05]", "(0.2,0.25]", "(0.4,0.45]", "(0.6,0.65]", "(0.8,0.85]"),
    labels = c("0", "0.2", "0.4", "0.6", "0.8")
  ) +
  theme_gray(base_size = 14)

ggsave("figures/qcrsf-heatmap.png", device = ragg::agg_png, width = 8, height = 6)

tmp <- qcs3_hcpya |>
  select(-window_start) |>
  filter(iter > 0) |>
  group_nest(sub, max_fd, filtered, cleaned) |>
  left_join(ref) |>
  mutate(q = purrr::map2_dbl(data, r, ~ ecdf(.x$r)(.y))) |>
  select(-data) |>
  mutate(max_fd = cut(max_fd, breaks = c(seq(0, 1, by = 0.05), Inf))) |>
  mutate(n_subs = n_distinct(sub), .by = c(max_fd)) |>
  filter(q < 0.05) |>
  distinct(sub, max_fd, filtered, cleaned) |>
  count(max_fd, filtered, cleaned) |>
  mutate(cdf = n / n_subs)

tmp |>
  ggplot(aes(x = max_fd, y = cdf, color = filtered, group = interaction(filtered, cleaned))) +
  geom_line(aes(linetype = cleaned))
