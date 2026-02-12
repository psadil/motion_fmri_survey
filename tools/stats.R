library(dplyr)

# there were XXXXXXX scans with motion data
targets::tar_load(by_run)

nrow(by_run)

# Across all datasets, a point higher corresponds to about XXXXX higher average motion.

targets::tar_load(c(by_run, demographics))
all_by_run <- by_run |>
  filter(!filtered) |>
  left_join(distinct(demographics, dataset, sub, ses, age, sex_gender, bmi))

cor.test(
  all_by_run$bmi,
  all_by_run$loc,
  use = "complete.obs",
  method = "spear"
) |>
  report::report()

# motion across sex

t.test(loc ~ sex_gender, data = all_by_run) |> report::report()

# In a linear model predicting average FD by sex, age, and their interaction

lm(loc ~ sex_gender * age, data = all_by_run) |> report::report()

# compare lost filtered
targets::tar_load(datasets)
source("R/utils.R")

n_subs <- by_run |>
  dplyr::semi_join(dplyr::distinct(datasets, dataset)) |>
  dplyr::count(dataset, task, ses, scan, filtered, name = "n_sub")

by_prop <- .summarize_by_prop(dataset, threshold = 0.2) |>
  dplyr::filter(lost > prop_thresh) |>
  dplyr::distinct(dataset, task, ses, sub, scan, filtered)
by_prop2 <- by_prop |>
  dplyr::count(dataset, task, ses, scan, filtered) |>
  dplyr::right_join(n_subs) |>
  dplyr::mutate(n = dplyr::if_else(is.na(n), 0, n)) |>
  dplyr::mutate(prop = n / n_sub) |>
  dplyr::select(-n, -n_sub)

by_max <- .summarize_by_max(dataset) |>
  dplyr::filter(max_fd > max_fd_thresh) |>
  dplyr::distinct(dataset, task, ses, sub, scan, filtered)
by_max2 <- by_max |>
  dplyr::count(dataset, task, ses, scan, filtered) |>
  dplyr::right_join(n_subs) |>
  dplyr::mutate(n = dplyr::if_else(is.na(n), 0, n)) |>
  dplyr::mutate(max = n / n_sub) |>
  dplyr::select(-n, -n_sub)

by_avg <- by_run |>
  dplyr::filter(loc > mfd_thresh) |>
  dplyr::distinct(dataset, task, ses, sub, scan, filtered)
by_avg2 <- by_avg |>
  dplyr::count(dataset, task, ses, scan, filtered) |>
  dplyr::right_join(n_subs) |>
  dplyr::mutate(n = dplyr::if_else(is.na(n), 0, n)) |>
  dplyr::mutate(avg = n / n_sub) |>
  dplyr::select(-n, -n_sub)

dplyr::bind_rows(by_max, by_avg, by_prop) |>
  dplyr::distinct() |>
  dplyr::count(dataset, task, ses, scan, filtered) |>
  dplyr::right_join(
    n_subs,
    by = dplyr::join_by(dataset, task, ses, scan, filtered)
  ) |>
  dplyr::mutate(n = dplyr::if_else(is.na(n), 0, n)) |>
  dplyr::mutate(
    lost = n / n_sub,
    lower = qbeta(0.025, 1 / 2 + n, n_sub - n + 1 / 2),
    upper = qbeta(0.975, 1 / 2 + n, n_sub - n + 1 / 2)
  ) |>
  dplyr::left_join(
    by_prop2,
    by = dplyr::join_by(dataset, task, ses, scan, filtered)
  ) |>
  dplyr::left_join(
    by_max2,
    by = dplyr::join_by(dataset, task, ses, scan, filtered)
  ) |>
  dplyr::left_join(
    by_avg2,
    by = dplyr::join_by(dataset, task, ses, scan, filtered)
  )


glm(
  lost ~ filtered * dataset,
  family = "binomial",
  data = filter(lost, type == "strict")
) |>
  report::report()
