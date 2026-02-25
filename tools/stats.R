library(dplyr)
library(tidyr)

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

# 1. Demonstrate downstream bias in brain–behavior associations

targets::tar_load(c(lost_strict2))

ukb_demo <- targets::tar_read(demographics) |>
  filter(dataset == "ukb", forcats::fct_match(ses, "2")) |>
  mutate(sub = as.double(sub))


bulk <- duckplyr::read_parquet_duckdb("data/ukb677207_bulk.parquet") |>
  select(sub = eid, etiv = f.26521.2.0) |>
  right_join(ukb_demo) |>
  collect()

cor.test(bulk$etiv, bulk$bmi) |> broom::tidy()


#' 2. Threshold sensitivity curves
#' done

# 3. Compute ICC

# across sessions

consistency <- by_run |>
  select(-sem) |>
  filter(stringr::str_detect(task, "rest"), scan == 1) |>
  group_nest(dataset, filtered) |>
  mutate(
    fit = purrr::map(
      data,
      ~ .x |>
        pivot_wider(names_from = ses, values_from = loc) |>
        select(-task, -scan, -sub) |>
        as.matrix() |>
        irr::icc(model = "t", type = "c")
    ),
    .lower = purrr::map_dbl(fit, purrr::pluck, "lbound"),
    .upper = purrr::map_dbl(fit, purrr::pluck, "ubound"),
    .estimate = purrr::map_dbl(fit, purrr::pluck, "value")
  )

# NOTE: updated from irr::icc to psych::ICC, but haven't updated manuscript
consistency <- by_run |>
  select(-sem) |>
  filter(
    stringr::str_detect(task, "rest"),
    scan == 1,
    dataset %in% c("ukb", "abcd")
  ) |>
  group_nest(dataset, filtered) |>
  mutate(
    fit = purrr::map(
      data,
      ~ .x |>
        pivot_wider(names_from = ses, values_from = loc) |>
        select(-task, -scan, -sub) |>
        as.matrix() |>
        psych::ICC() |>
        purrr::pluck("results") |>
        as_tibble()
    )
  ) |>
  select(-data) |>
  unnest(fit)


# discriminability by session
# dd <- by_run |>
#   select(-sem) |>
#   unite(scan, task, scan) |>
#   group_nest(dataset, filtered) |>
#   filter(dataset %in% c("abcd", "ukb")) |>
#   mutate(
#     data = purrr::map(
#       data,
#       ~ .x |>
#         pivot_wider(names_from = scan, values_from = loc) |>
#         na.omit()
#     ),
#     dis = purrr::map(
#       data,
#       ~ mgc::discr.stat(
#         select(.x, -ses, -sub),
#         .x$sub,
#         no_cores = 4
#       ),
#       .progress = TRUE
#     )
#   )

# icc by scan (rest)
# dd <- by_run |>
#   select(-sem) |>
#   filter(!(dataset %in% c("ukb", "spacetop"))) |>
#   filter(task == "rest") |>
#   group_nest(dataset, filtered, ses) |>
#   mutate(
#     fit = purrr::map(
#       data,
#       ~ .x |>
#         pivot_wider(names_from = scan, values_from = loc) |>
#         na.omit() |>
#         select(`1`:`4`) |>
#         as.matrix() |>
#         irr::icc("t")
#     ),
#     .lower = purrr::map_dbl(fit, purrr::pluck, "lbound"),
#     .upper = purrr::map_dbl(fit, purrr::pluck, "ubound"),
#     .estimate = purrr::map_dbl(fit, purrr::pluck, "value")
#   )

dd <- by_run |>
  select(-sem) |>
  filter(!(dataset %in% c("ukb", "spacetop"))) |>
  filter(task == "rest") |>
  group_nest(dataset, filtered, ses) |>
  mutate(
    fit = purrr::map(
      data,
      ~ .x |>
        pivot_wider(names_from = scan, values_from = loc) |>
        na.omit() |>
        select(`1`:`4`) |>
        as.matrix() |>
        psych::ICC() |>
        purrr::pluck("results") |>
        as_tibble()
    )
  ) |>
  select(-data) |>
  unnest(fit)


# 4. Logistic model of exclusion probability

targets::tar_load(c(lost_strict2, demographics))

d <- lost_strict2 |>
  filter(
    !filtered,
    task == "rest",
    scan == 1,
    ses %in% c("baseline", "2", "1")
  ) |>
  left_join(distinct(
    demographics,
    dataset,
    sub,
    ses,
    age,
    sex = sex_gender,
    bmi
  )) |>
  na.omit()

fit <- lme4::glmer(
  exclude ~ poly(age, 2) + sex * bmi + (1 | ses:dataset) + (1 | sub),
  data = d,
  family = "binomial",
  control = lme4::glmerControl(autoscale = TRUE)
)

fit <- gamm4::gamm4(
  exclude ~ s(age) + sex * bmi + dataset,
  random = ~ (1 | sub),
  data = d,
  family = "binomial"
)


# fit <- glm(
#   exclude ~ poly(age, 2) + sex + bmi + dataset,
#   data = d,
#   family = "binomial"
# )

fit <- glmmTMB::glmmTMB(
  exclude ~ poly(age, 3) + sex * bmi + (1 | ses:dataset) + (1 | sub),
  data = d,
  family = "binomial",
)

fit <- glmmTMB::glmmTMB(
  exclude ~ s(age) + sex * bmi + (1 | ses:dataset) + (1 | sub),
  data = d,
  family = "binomial",
  REML = TRUE,
)

# e <- emmeans::emmeans(fit, ~ poly(age, 2) * sex + bmi + dataset)

broom.mixed::tidy(fit)
broom.mixed::tidy(fit, exponentiate = TRUE)

# library(brms)

# f <- bf(exclude ~ s(age) + sex + bmi + (1 | dataset) + (1 | sub))
# bp <- c(prior(cauchy(0, 2.5), class = b))

# fit_b <- brm(
#   formula = f,
#   data = filter(d, stringr::str_detect(dataset, "hcp") | dataset == "ukb"),
#   family = bernoulli(),
#   cores = 4,
#   chains = 4,
#   save_model = "tmp/model.stan",
#   prior = bp
# )

# power analysis
targets::tar_load(c(lost, by_run))
library(ggplot2)

by_run |>
  filter(scan == 1, task == "rest") |>
  count(dataset, ses, task, scan, filtered) |>
  left_join(lost) |>
  mutate(nn = n * (1 - lost)) |>
  rename(full = n, Filtered = filtered) |>
  select(-lower, -upper, -prop, -max, -avg, -lost) |>
  pivot_wider(names_from = type, values_from = nn) |>
  pivot_longer(full:strict, values_to = "n") |>
  mutate(
    p = pwr::pwr.r.test(n = n, r = 0.1)$power,
    dataset = stringr::str_to_upper(dataset),
    ses = stringr::str_replace(ses, "^[[:digit:]]+", " ")
  ) |>
  ggplot(aes(x = name, y = p)) +
  geom_col(aes(fill = Filtered), position = "dodge") +
  facet_wrap(~ dataset + ses) +
  scale_fill_viridis_d(option = "turbo") +
  ylab("Statistical Power")
