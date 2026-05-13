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


# UKB subset info:

targets::tar_load(c(by_run, ukb_subs))

by_run |>
  filter(dataset == "ukb", !filtered) |>
  mutate(sub = as.numeric(sub)) |>
  filter(sub %in% ukb_subs, task == "rest", ses == "2") |>
  summarise(avg = mean(loc), s = sd(loc))


# demographic table of high vs. low movers
library(dplyr)
library(tidyr)
library(table1)
library(ggplot2)

targets::tar_load(c(by_run, demographics, lost_strict2))

splits <- by_run |>
  filter(
    task == "rest",
    !filtered,
    scan == 1,
    dataset %in% c("abcd", "hcpya"),
    ses %in% c("Baseline", 1)
  ) |>
  mutate(group = if_else(loc > median(loc), "high", "low"), .by = c(dataset)) |>
  select(dataset, sub, group, loc)

toshow <- demographics |>
  inner_join(splits) |>
  filter(ses %in% c(1, "Baseline")) |>
  rename(sex = sex_gender) |>
  mutate(sex = forcats::fct_drop(sex))

table1(~ age + bmi + sex | dataset + group, toshow)

library(gtsummary)

toshow |>
  select(dataset, Age = age, BMI = bmi, Sex = sex, Group = group) |>
  na.omit() |>
  gtsummary::tbl_strata(
    strata = dataset,
    ~ .x |>
      tbl_summary(by = Group) |>
      modify_header(all_stat_cols() ~ "**{level}**")
  ) |>
  as_gt() |>
  gt::gtsave("figures/highlowdemo.tex")


toshow2 <- readr::read_csv(
  "data/human-connectome-project-restricted/RESTRICTED_martin_2_5_2024_10_18_28.csv"
) |>
  select(sub = Subject, starts_with("ASR_"), starts_with("DSM_")) |>
  select(!ends_with("Pct") & !ends_with("Sum")) |>
  mutate(sub = as.character(sub)) |>
  rename_with(
    ~ .x |>
      stringr::str_remove("ASR_") |>
      stringr::str_remove("DSM_") |>
      stringr::str_remove("_T")
  ) |>
  inner_join(splits) |>
  select(-dataset)

table1(
  ~ Witd +
    Soma +
    Thot +
    Attn +
    Aggr +
    Rule +
    Intr +
    Intn +
    Extn +
    Totp +
    Depr +
    Anxi +
    Somp +
    Avoid +
    Adh +
    Antis |
    group,
  toshow2
)


toshow2 |>
  select(-ends_with("Raw")) |>
  pivot_longer(c(-sub, -group, -loc)) |>
  group_nest(name) |>
  mutate(
    data = purrr::map(
      data,
      ~ lm(value ~ loc, .x) |> broom::tidy(conf.int = TRUE)
    )
  ) |>
  unnest(data) |>
  filter(stringr::str_detect(term, "Inter", TRUE)) |>
  arrange(estimate) |>
  mutate(name = factor(name, levels = .data$name)) |>
  ggplot(aes(y = name)) +
  geom_point(aes(x = estimate)) +
  geom_segment(aes(x = conf.low, xend = conf.high)) +
  ylab("Psychiatric and Life Function Variable") +
  xlab("Estimate")

ggsave("figures/hcp-psych.png", width = 4, height = 4, device = ragg::agg_png)

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


library(gtsummary)
targets::tar_load(c(demographics))

plot_demo_tbl <- function(by_run, demographics, ds) {
  by_run |>
    filter(dataset == ds) |>
    distinct(sub, ses, dataset) |>
    left_join(demographics) |>
    select(sub, ses, Age = age, BMI = bmi, Sex = sex_gender) |>
    mutate(ses = forcats::fct_drop(ses), Sex = forcats::fct_drop(Sex)) |>
    tbl_summary(by = ses, include = c(-sub)) |>
    as_gt() |>
    gt::gtsave(glue::glue("figures/demographics-{ds}.tex"))
}

purrr::walk(unique(by_run$dataset), ~ plot_demo_tbl(by_run, demographics, .x))

# spectrum peaks table
get_limits <- function(d, lower_exclude = 0.15, upper_exclude = 0.6) {
  d |>
    filter(between(freq, lower_exclude, upper_exclude), param == "trans_y") |>
    slice_max(order_by = pxx, n = 1, by = c(sub, ses, ped, task)) |>
    summarise(freq = median(freq), .by = c(sub, param, task)) |>
    summarise(
      lower = median(freq) |> round(2),
      upper = quantile(freq, 0.75) |> round(2)
    )
}

targets::tar_load(c(
  demographics,
  hcpa_spectrum,
  hcpd_spectrum,
  hcpya_spectrum,
  spacetop_spectrum
))

get_limits(hcpa_spectrum)
hcpd_spectrum |>
  left_join(
    demographics |> filter(dataset == "hcpd") |> distinct(sub, ses, age)
  ) |>
  filter(age < 8) |>
  get_limits()

hcpd_spectrum |>
  left_join(
    demographics |> filter(dataset == "hcpd") |> distinct(sub, ses, age)
  ) |>
  filter(age >= 8) |>
  get_limits()

get_limits(hcpya_spectrum)
get_limits(spacetop_spectrum)


abcd <- arrow::read_parquet(
  "/Users/psadil/git/manuscripts/motion/data/motion/derivatives/abcd_spectrum.parquet"
) |>
  filter(between(freq, 0.2, 0.6), param == "trans_y") |>
  slice_max(order_by = pxx, n = 1, by = c(sub, ses, task)) |>
  summarise(
    lower = median(freq),
    upper = quantile(freq, 0.75),
    .by = c(task, ses)
  )

ukb <- arrow::read_parquet(
  "/Users/psadil/git/manuscripts/motion/data/motion/derivatives/ukb_spectrum.parquet"
) |>
  filter(between(freq, 0.2, 0.6), param == "trans_y") |>
  slice_max(order_by = pxx, n = 1, by = c(sub, ses, task)) |>
  summarise(freq = median(freq), .by = c(sub, task)) |>
  summarise(lower = median(freq), upper = quantile(freq, 0.75), .by = c(task))


# comparison of ABCD with HCPD

library(dplyr)
library(tidyr)
targets::tar_load(c(by_run, demographics))


bpm <- readr::read_csv("data/tabular/core/mental-health/mh_t_bpm.csv") |>
  select(
    sub = src_subject_id,
    ses = eventname,
    bpm_t_scr_external_t,
    bpm_t_scr_internal_t,
    bpm_t_scr_attention_t
  ) |>
  filter(
    !(ses %in% c("1_year_follow_up_y_arm_1", "3_year_follow_up_y_arm_1"))
  ) |>
  mutate(sub = stringr::str_remove(sub, "_")) |>
  convert_abcd_ses() |>
  pivot_longer(bpm_t_scr_external_t:bpm_t_scr_attention_t) |>
  na.omit() |>
  summarise(value = mean(value), .by = c(sub, name)) |>
  pivot_wider()


d <- by_run |>
  filter(
    dataset %in% c("hcpd", "abcd"),
    !filtered,
    stringr::str_detect(task, "rest")
  ) |>
  left_join(distinct(demographics, dataset, sub, age, ses, bmi, sex)) |>
  left_join(bpm) |>
  mutate(across(contains("bpm_t_scr"), ~ if_else(is.na(.x), 0, .x))) |>
  mutate(
    external = bpm_t_scr_external_t > 65,
    internal = bpm_t_scr_internal_t > 65,
    attention = bpm_t_scr_attention_t > 65,
  )

fit <- lme4::lmer(loc ~ dataset + age * sex_gender * bmi + (1 | sub), data = d)

fit |> report::report()

d2 <- by_run |>
  filter(dataset %in% c("hcpd", "abcd"), !filtered) |>
  left_join(distinct(demographics, dataset, sub, age, ses, bmi, sex_gender)) |>
  filter(age < 10)

fit2 <- lme4::lmer(
  loc ~ dataset + age * sex_gender * bmi + (1 | sub),
  data = d2
)

fit2 |> report::report()

fit3 <- lme4::lmer(
  loc ~ dataset + age * bmi * sex_gender + (1 | sub),
  data = d |>
    mutate(
      young = age < 8,
      dataset = interaction(dataset, ses, young, drop = TRUE)
    )
)

fit3 |> report::report()


fit4 <- lme4::lmer(
  loc ~ dataset + dataset:age + age * bmi * sex + (1 | sub),
  data = d
)

fit4 |> report::report()

fit4 |>
  report::report() |>
  as.data.frame() |>
  insight::format_table() |>
  insight::export_table()


fit5 <- lme4::lmer(
  loc ~ external +
    internal +
    attention +
    dataset +
    dataset:age +
    age * bmi * sex +
    (1 | sub),
  data = d
)

fit5 |> report::report()


# abcd lost by race/ethnicity
dd <- readr::read_csv(
  "data/tabular/core/abcd-general/abcd_p_demo.csv",
  show_col_types = FALSE
) |>
  select(sub = src_subject_id, ses = eventname, race_ethnicity) |>
  filter((ses == "baseline_year_1_arm_1")) |>
  mutate(sub = stringr::str_remove(sub, "_")) |>
  convert_abcd_ses() |>
  mutate(
    race_ethnicity = dplyr::recode_values(
      race_ethnicity,
      1 ~ "White",
      2 ~ "Black",
      3 ~ "Hispanic",
      4 ~ "Asian",
      5 ~ "Other"
    )
  )


d <- targets::tar_read(lost_strict2) |>
  filter(
    !filtered,
    ses == "Baseline",
    dataset == "abcd",
    task == "rest",
    scan == 1
  ) |>
  select(sub, exclude) |>
  left_join(dd)

d |> gtsummary::tbl_summary(include = c(race_ethnicity))

d |> gtsummary::tbl_summary(include = c(race_ethnicity), by = exclude)
