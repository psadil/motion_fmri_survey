library(dplyr)
library(ggplot2)
library(tidyr)
library(ggdist)
library(patchwork)

targets::tar_load(c(ukb_reg, hcp))

add_fd <- function(.data, ..., radius = 50) {
  .data |>
    group_by(...) |>
    arrange(t, .by_group = TRUE) |>
    mutate(
      across(matches("[trans|rot]_[xyz]"), \(x) abs(lag(x) - x))
    ) |>
    ungroup() |>
    na.omit() |>
    dplyr::mutate(
      dplyr::across(tidyselect::starts_with("rot"), \(x) x * radius),
      fd = rot_x + rot_y + rot_z + trans_x + trans_y + trans_z
    )
}
fd_max_y <- 0.4

hcp_tmp <- hcp |>
  select(-ends_with("derivative1")) |>
  add_fd(sub, task, ped) |>
  summarise(
    mean = mean(fd),
    .by = c(task, ped, ses, sub)
  ) |>
  mutate(ses = glue::glue("Ses {ses}"))

ukb_tmp <- ukb_reg |>
  add_fd(sub, ses) |>
  summarise(
    mean = mean(fd),
    .by = c(ses, sub)
  ) |>
  mutate(ses = glue::glue("Ses {ses}"))

a <- hcp_tmp |>
  ggplot(aes(x = mean, color = ped, fill = ped)) +
  facet_wrap(~task) +
  stat_slab(aes(thickness = after_stat(pdf * n)), scale = 0.7, alpha = 0.5) +
  stat_dotsinterval(
    side = "bottom",
    scale = 0.7,
    slab_linewidth = NA,
    quantiles = 100,
    alpha = 0.5
  ) +
  ylab("Repetition Time (grouped, seconds)") +
  xlab("Framewise Displacement Mean (mm)")

b <- ukb_tmp |>
  ggplot(aes(x = mean)) +
  stat_slab(aes(thickness = after_stat(pdf * n)), scale = 0.7, alpha = 0.5) +
  stat_dotsinterval(
    side = "bottom",
    scale = 0.7,
    slab_linewidth = NA,
    quantiles = 100,
    alpha = 0.5
  ) +
  ylab("Repetition Time (grouped, seconds)") +
  xlab("Framewise Displacement Mean (mm)")
