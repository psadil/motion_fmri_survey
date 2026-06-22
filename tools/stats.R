library(dplyr)
library(tidyr)


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


# 4. Logistic model of exclusion probability

targets::tar_load(c(lost_strict2, demographics))

d <- lost_strict2 |>
  filter(
    !filtered,
    task == "rest",
    scan == 1,
    ses %in% c("Baseline", "2", "1")
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

# UKB subset info:

targets::tar_load(c(by_run, ukb_subs))

by_run |>
  filter(dataset == "ukb", !filtered) |>
  mutate(sub = as.numeric(sub)) |>
  filter(sub %in% ukb_subs, task == "rest", ses == "2") |>
  summarise(avg = mean(loc), s = sd(loc))


nih <- readr::read_csv("data/tabular/core/neurocognition/nc_y_nihtb.csv") |>
  select(sub = src_subject_id, ses = eventname, ends_with("agecorrected")) |>
  rename_with(
    ~ stringr::str_remove(.x, "nihtbx_") |> stringr::str_remove("_agecorrected")
  ) |>
  mutate(
    ses = replace_values(
      ses,
      "baseline_year_1_arm_1" ~ "Baseline",
      "2_year_follow_up_y_arm_1" ~ "Year2",
      "4_year_follow_up_y_arm_1" ~ "Year4"
    ),
    sub = stringr::str_remove(sub, "_")
  )

xsx <- readr::read_csv("data/tabular/core/imaging/mri_y_rsfmr_cor_gp_gp.csv") |>
  rename(sub = src_subject_id, ses = eventname) |>
  mutate(
    ses = replace_values(
      ses,
      "baseline_year_1_arm_1" ~ "Baseline",
      "2_year_follow_up_y_arm_1" ~ "Year2",
      "4_year_follow_up_y_arm_1" ~ "Year4"
    ),
    sub = stringr::str_remove(sub, "_")
  )

d_s <- left_join(nih, xsx) |>
  inner_join(
    filter(
      lost_strict2,
      dataset == "abcd",
      task == "rest",
      scan == "1",
      !filtered
    ) |>
      distinct()
  )

predictor_vars <- nih |> select(-sub, -ses) |> names()
outcome_vars <- xsx |> select(-sub, -ses) |> names()

get_cors <- function(d, predictor_vars, outcome_vars) {
  cors <- expand_grid(predictor = predictor_vars, outcome = outcome_vars) |>
    mutate(
      # Use map2_dbl to iterate over the two columns simultaneously
      correlation = purrr::map2(
        .x = predictor,
        .y = outcome,
        ~ cor.test(d[[.x]], d[[.y]], method = "spear") |>
          broom::tidy() |>
          select(-method, -alternative),
        .progress = TRUE
      )
    ) |>
    unnest(correlation)
}

nolost <- get_cors(d_s, predictor_vars, outcome_vars)
strict <- get_cors(filter(d_s, exclude), predictor_vars, outcome_vars)

bind_rows(list("none" = nolost, "strict" = strict), .id = "group") |>
  select(-statistic, -p.value) |>
  pivot_wider(names_from = group, values_from = estimate) |>
  ggplot(aes(x = none, y = strict)) +
  geom_point() +
  geom_abline() +
  geom_smooth(method = "lm")

ggsave("figures/strict-vs-none-bwas.png", device = ragg::agg_png)
