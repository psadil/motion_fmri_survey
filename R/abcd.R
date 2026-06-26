get_abcd_entities <- function(abcd_source) {
  tibble::tibble(src = abcd_source) |>
    dplyr::mutate(
      sub = stringr::str_extract(src, "(?<=sub-)([[:digit:]]|[[:alpha:]])+"),
      ses = stringr::str_extract(src, "(?<=ses-)([[:digit:]]|[[:alpha:]])+"),
      task = stringr::str_extract(src, "(?<=task-)[[:alpha:]]+"),
      run = stringr::str_extract(src, "(?<=run-)[[:digit:]]+")
    ) |>
    do_casting() |>
    convert_abcd_ses()
}

get_abcd <- function(abcd_source, abcd_exclusion, abcd_phenotypes) {
  entities <- get_abcd_entities(abcd_source)

  duckplyr::read_file_duckdb(
    abcd_source,
    prudence = "lavish",
    table_function = "read_csv_auto",
    options = list(filename = "src")
  ) |>
    dplyr::group_nest(src) |>
    dplyr::mutate(
      data = purrr::map(data, ~ dplyr::mutate(.x, t = dplyr::row_number()))
    ) |>
    tidyr::unnest(data) |>
    dplyr::mutate(t = t - 1) |>
    dplyr::filter(t > 0) |>
    dplyr::left_join(entities, by = dplyr::join_by(src)) |>
    dplyr::select(-src) |>
    dplyr::anti_join(abcd_exclusion, dplyr::join_by(sub, ses)) |>
    dplyr::mutate(time = t * 0.8, ped = "AP", scan = run) |>
    dplyr::collect() |>
    exclude_bad_abcd_scan(abcd_source, abcd_phenotypes) |>
    truncate_to_modal_lengths()
}

convert_abcd_ses <- function(.d) {
  .d |>
    dplyr::mutate(
      ses = dplyr::replace_values(
        ses,
        "baseline_year_1_arm_1" ~ "Baseline",
        "2_year_follow_up_y_arm_1" ~ "Year2",
        "4_year_follow_up_y_arm_1" ~ "Year4",
        "ses-00A" ~ "Baseline",
        "ses-02A" ~ "Year2",
        "ses-04A" ~ "Year4",
        "ses-06A" ~ "Year6",
        "00A" ~ "Baseline",
        "02A" ~ "Year2",
        "04A" ~ "Year4",
        "06A" ~ "Year6",
        "baselineYear1Arm1" ~ "Baseline",
        "2YearFollowUpYArm1" ~ "Year2",
        "4YearFollowUpYArm1" ~ "Year4"
      )
    )
}

get_abcd_events <- function(abcd_design) {
  events <- duckplyr::read_parquet_duckdb(abcd_design, prudence = "lavish")

  cue_design <- events |>
    dplyr::filter(task == "nback") |>
    dplyr::filter(trial_type == "cue") |>
    dplyr::arrange(onset) |>
    dplyr::collect() |>
    dplyr::mutate(event = 1:dplyr::n(), .by = c(run, task, sub, ses)) |>
    dplyr::summarise(
      onset = mean(onset),
      duration = mean(duration),
      .by = c(event, task)
    ) |>
    dplyr::mutate(type = "cue")

  nback_design <- events |>
    dplyr::filter(stringr::str_detect(task, "nback")) |>
    dplyr::filter(stringr::str_detect(trial_type, "back")) |>
    dplyr::arrange(onset) |>
    dplyr::collect() |>
    dplyr::mutate(event = 1:dplyr::n(), .by = c(run, task, sub, ses)) |>
    dplyr::summarise(
      onset = mean(onset),
      duration = mean(duration),
      .by = c(event, task)
    ) |>
    dplyr::mutate(type = "nback")

  nback <- dplyr::bind_rows(cue_design, nback_design)

  # SST timing varies by performance

  # sst <- events |>
  #   dplyr::filter(task == "sst") |>
  #   dplyr::filter(stringr::str_detect(trial_type, "go")) |>
  #   dplyr::arrange(onset) |>
  #   dplyr::collect() |>
  #   dplyr::mutate(
  #     event = 1:dplyr::n(),
  #     .by = c(run, task, sub, ses)
  #   ) |>
  #   dplyr::summarise(
  #     onset = mean(onset),
  #     duration = mean(duration),
  #     .by = c(event, task)
  #   ) |>
  #   dplyr::mutate(type = "go/nogo")

  mid_antic <- events |>
    dplyr::filter(stringr::str_detect(trial_type, "antic")) |>
    dplyr::arrange(onset) |>
    dplyr::collect() |>
    dplyr::mutate(event = 1:dplyr::n(), .by = c(run, task, sub, ses)) |>
    dplyr::summarise(
      onset = mean(onset),
      duration = mean(duration),
      .by = c(event, task)
    ) |>
    dplyr::mutate(type = "antic")

  mid_fb <- events |>
    dplyr::filter(task == "mid") |>
    dplyr::filter(stringr::str_detect(trial_type, "feedback")) |>
    dplyr::arrange(onset) |>
    dplyr::collect() |>
    dplyr::mutate(event = 1:dplyr::n(), .by = c(run, task, sub, ses)) |>
    dplyr::summarise(
      onset = mean(onset),
      duration = mean(duration),
      .by = c(event, task)
    ) |>
    dplyr::mutate(type = "feedback")

  mid <- dplyr::bind_rows(mid_antic, mid_fb)

  # dplyr::bind_rows(nback, sst, mid)
  dplyr::bind_rows(nback, mid)
}

read_abcd_phenotype <- function(src) {
  d <- duckplyr::read_parquet_duckdb(src, prudence = "lavish") |>
    dplyr::rename(sub = participant_id) |>
    dplyr::mutate(sub = stringr::str_remove(sub, "sub-"))

  if ("session_id" %in% names(d)) {
    d <- d |>
      dplyr::rename(ses = "session_id") |>
      dplyr::filter(ses %in% c("ses-00A", "ses-02A", "ses-04A", "ses-06A")) |>
      convert_abcd_ses()
  }
  d
}

get_abcd_demographics <- function(abcd_source) {
  abcd_subs <- get_abcd_subs_from_srcs(abcd_source)

  serial <- read_abcd_phenotype(
    "data/abcc-4-0-0/rawdata/phenotype/mr_y_adm__info.parquet"
  ) |>
    dplyr::select(
      sub,
      ses,
      mr_y_adm__info_dt,
      deviceserialnumber = mr_y_adm__info__dev_serial,
      manufacturer = mr_y_adm__info__dev_manufact
    ) |>
    dplyr::mutate(
      manufacturer = dplyr::case_when(
        manufacturer == 1 ~ "GE",
        manufacturer == 2 ~ "Philips",
        manufacturer == 3 ~ "Siemens",
        TRUE ~ NA_character_
      )
    )

  age <- read_abcd_phenotype(
    "data/abcc-4-0-0/rawdata/phenotype/ab_p_demo.parquet"
  ) |>
    dplyr::select(sub, ses, age = ab_p_demo_age)

  sex <- read_abcd_phenotype(
    "data/abcc-4-0-0/rawdata/phenotype/ab_g_stc.parquet"
  ) |>
    dplyr::select(sub, sex = ab_g_stc__cohort_sex) |>
    dplyr::mutate(
      sex = dplyr::recode_values(
        sex,
        "1" ~ "Male",
        "2" ~ "Female",
        default = NA_character_
      )
    )

  bmi <- read_abcd_phenotype(
    "data/abcc-4-0-0/rawdata/phenotype/ph_y_anthr.parquet"
  ) |>
    dplyr::select(
      sub,
      ses,
      weight = ph_y_anthr__weight_mean,
      height = ph_y_anthr__height_mean
    ) |>
    dplyr::mutate(
      weight = weight * 0.4535924,
      height = height * 0.0254,
      bmi = weight / height^2
    ) |>
    dplyr::select(sub, ses, bmi)

  dplyr::full_join(serial, bmi, by = dplyr::join_by(sub, ses)) |>
    dplyr::full_join(sex, by = dplyr::join_by(sub)) |>
    dplyr::full_join(age, by = dplyr::join_by(sub, ses)) |>
    dplyr::filter(ses %in% c("Baseline", "Year2", "Year4", "Year6")) |>
    dplyr::filter(sub %in% abcd_subs) |>
    dplyr::select(
      sub,
      ses,
      age,
      sex,
      deviceserialnumber,
      bmi,
      manufacturer,
      mr_y_adm__info_dt
    ) |>
    dplyr::collect()
}

get_abcd_exclusion_crit <- function(src, rule) {
  read_abcd_phenotype(src) |> dplyr::filter({{ rule }})
}


get_abcd_excl_rest <- function(srcs) {
  # https://docs.abcdstudy.org/latest/documentation/imaging/type_qc.html#rs-fmri-data-recommended-for-inclusion
  # everything except censoring

  tfMRI_series_passed_rawQC <- get_abcd_exclusion_crit(
    fs::path(srcs, "mr_y_qc__raw__rsfmri.parquet"),
    !(mr_y_qc__raw__rsfmri__pass__qc__comp_count > 0)
  )

  T1_series_passed_rawQC <- get_abcd_exclusion_crit(
    fs::path(srcs, "mr_y_qc__raw__smri__t1.parquet"),
    !(mr_y_qc__raw__smri__t1__pass__qc__comp_count > 0)
  )

  fMRI_B0_unwarp_available <- get_abcd_exclusion_crit(
    fs::path(srcs, "mr_y_qc__post__aut.parquet"),
    !(mr_y_qc__post__aut__fmri__b0__unwarp_indicator == 1)
  )

  FreeSurfer_QC_not_failed <- get_abcd_exclusion_crit(
    fs::path(srcs, "mr_y_qc__post__man__fsurf.parquet"),
    (mr_y_qc__post__man__fsurf_score == 0)
  )

  fMRI_manual_post_processing_QC_not_failed <- get_abcd_exclusion_crit(
    fs::path(srcs, "mr_y_qc__post__man__fmri.parquet"),
    (mr_y_qc__post__man__fmri_score == 0)
  )

  fMRI_registration_to_T1w <- get_abcd_exclusion_crit(
    fs::path(srcs, "mr_y_qc__post__aut.parquet"),
    mr_y_qc__post__aut__fmri__rigid_score >= 19
  )

  fMRI_dorsal_cutoff_score <- get_abcd_exclusion_crit(
    fs::path(srcs, "mr_y_qc__post__aut.parquet"),
    mr_y_qc__post__aut__fmri__fov__cutoff__dorsal_max >= 65
  )

  fMRI_ventral_cutoff_score <- get_abcd_exclusion_crit(
    fs::path(srcs, "mr_y_qc__post__aut.parquet"),
    mr_y_qc__post__aut__fmri__fov__cutoff__ventral_max >= 60
  )

  Derived_results_exist <- get_abcd_exclusion_crit(
    fs::path(srcs, "mr_y_rsfmri__corr__gpnet.parquet"),
    is.na(mr_y_rsfmri__corr__gpnet__def__sal_mean)
  )

  purrr::reduce(
    .x = list(
      tfMRI_series_passed_rawQC,
      T1_series_passed_rawQC,
      fMRI_B0_unwarp_available,
      FreeSurfer_QC_not_failed,
      fMRI_manual_post_processing_QC_not_failed,
      fMRI_registration_to_T1w,
      fMRI_dorsal_cutoff_score,
      fMRI_ventral_cutoff_score,
      Derived_results_exist
    ),
    .f = dplyr::full_join,
    dplyr::join_by(sub, ses),
  ) |>
    dplyr::distinct(sub, ses) |>
    dplyr::mutate(task = "rest") |>
    dplyr::collect()
}


get_abcd_excl_mid <- function(srcs) {
  # https://docs.abcdstudy.org/latest/documentation/imaging/type_qc.html#mid-task-fmri-data-recommended-for-inclusion
  # everything except censoring

  tfMRI_series_passed_rawQC <- get_abcd_exclusion_crit(
    fs::path(srcs, "mr_y_qc__raw__tfmri__mid.parquet"),
    !(mr_y_qc__raw__tfmri__mid__pass__qc__comp_count > 0)
  )

  T1_series_passed_rawQC <- get_abcd_exclusion_crit(
    fs::path(srcs, "mr_y_qc__raw__smri__t1.parquet"),
    !(mr_y_qc__raw__smri__t1__pass__qc__comp_count > 0)
  )

  behavior_passed <- get_abcd_exclusion_crit(
    fs::path(srcs, "mr_y_tfmri__mid__beh.parquet"),
    !(mr_y_tfmri__mid__beh__qc_indicator == 1)
  )

  E_prime_timing_match_OR_ignore_E_prime_mismatch <- get_abcd_exclusion_crit(
    fs::path(srcs, "mr_y_qc__raw__tfmri__mid.parquet"),
    !(mr_y_qc__raw__tfmri__mid__eprime__match_indicator == 1 |
      mr_y_qc__raw__tfmri__mid__eprime__tdiff__ign_indicator < 1)
  )

  fMRI_B0_unwarp_available <- get_abcd_exclusion_crit(
    fs::path(srcs, "mr_y_qc__post__aut.parquet"),
    !(mr_y_qc__post__aut__fmri__b0__unwarp_indicator == 1)
  )

  FreeSurfer_QC_not_failed <- get_abcd_exclusion_crit(
    fs::path(srcs, "mr_y_qc__post__man__fsurf.parquet"),
    (mr_y_qc__post__man__fsurf_score == 0)
  )

  fMRI_manual_post_processing_QC_not_failed <- get_abcd_exclusion_crit(
    fs::path(srcs, "mr_y_qc__post__man__fmri.parquet"),
    (mr_y_qc__post__man__fmri_score == 0)
  )

  fMRI_registration_to_T1w <- get_abcd_exclusion_crit(
    fs::path(srcs, "mr_y_qc__post__aut.parquet"),
    mr_y_qc__post__aut__fmri__rigid_score >= 19
  )

  fMRI_dorsal_cutoff_score <- get_abcd_exclusion_crit(
    fs::path(srcs, "mr_y_qc__post__aut.parquet"),
    mr_y_qc__post__aut__fmri__fov__cutoff__dorsal_max >= 65
  )

  fMRI_ventral_cutoff_score <- get_abcd_exclusion_crit(
    fs::path(srcs, "mr_y_qc__post__aut.parquet"),
    mr_y_qc__post__aut__fmri__fov__cutoff__ventral_max >= 60
  )

  Derived_results_exist <- get_abcd_exclusion_crit(
    fs::path(srcs, "mr_y_tfmri__mid__arvn__aseg.parquet"),
    is.na(mr_y_tfmri__mid__arvn__aseg__cwm__lh_beta)
  )

  purrr::reduce(
    .x = list(
      tfMRI_series_passed_rawQC,
      T1_series_passed_rawQC,
      behavior_passed,
      E_prime_timing_match_OR_ignore_E_prime_mismatch,
      fMRI_B0_unwarp_available,
      FreeSurfer_QC_not_failed,
      fMRI_manual_post_processing_QC_not_failed,
      fMRI_registration_to_T1w,
      fMRI_dorsal_cutoff_score,
      fMRI_ventral_cutoff_score,
      Derived_results_exist
    ),
    .f = dplyr::full_join,
    dplyr::join_by(sub, ses),
  ) |>
    dplyr::distinct(sub, ses) |>
    dplyr::mutate(task = "mid") |>
    dplyr::collect()
}

get_abcd_excl_sst <- function(srcs) {
  # https://docs.abcdstudy.org/latest/documentation/imaging/type_qc.html#sst-task-fmri-data-recommended-for-inclusion
  # everything except censoring

  tfMRI_series_passed_rawQC <- get_abcd_exclusion_crit(
    fs::path(srcs, "mr_y_qc__raw__tfmri__sst.parquet"),
    !(mr_y_qc__raw__tfmri__sst__pass__qc__comp_count > 0)
  )

  T1_series_passed_rawQC <- get_abcd_exclusion_crit(
    fs::path(srcs, "mr_y_qc__raw__smri__t1.parquet"),
    !(mr_y_qc__raw__smri__t1__pass__qc__comp_count > 0)
  )

  behavior_passed <- get_abcd_exclusion_crit(
    fs::path(srcs, "mr_y_tfmri__sst__beh.parquet"),
    !(mr_y_tfmri__sst__beh__qc_indicator == 1)
  )

  task_had_no_glitch <- get_abcd_exclusion_crit(
    fs::path(srcs, "mr_y_tfmri__sst__beh.parquet"),
    !(mr_y_tfmri__sst__beh__coderr_indicator == 0)
  )

  E_prime_timing_match_OR_ignore_E_prime_mismatch <- get_abcd_exclusion_crit(
    fs::path(srcs, "mr_y_qc__raw__tfmri__sst.parquet"),
    !(mr_y_qc__raw__tfmri__sst__eprime__match_indicator == 1 |
      mr_y_qc__raw__tfmri__sst__eprime__tdiff__ign_indicator < 1)
  )

  fMRI_B0_unwarp_available <- get_abcd_exclusion_crit(
    fs::path(srcs, "mr_y_qc__post__aut.parquet"),
    !(mr_y_qc__post__aut__fmri__b0__unwarp_indicator == 1)
  )

  FreeSurfer_QC_not_failed <- get_abcd_exclusion_crit(
    fs::path(srcs, "mr_y_qc__post__man__fsurf.parquet"),
    (mr_y_qc__post__man__fsurf_score == 0)
  )

  fMRI_manual_post_processing_QC_not_failed <- get_abcd_exclusion_crit(
    fs::path(srcs, "mr_y_qc__post__man__fmri.parquet"),
    (mr_y_qc__post__man__fmri_score == 0)
  )

  fMRI_registration_to_T1w <- get_abcd_exclusion_crit(
    fs::path(srcs, "mr_y_qc__post__aut.parquet"),
    mr_y_qc__post__aut__fmri__rigid_score >= 19
  )

  fMRI_dorsal_cutoff_score <- get_abcd_exclusion_crit(
    fs::path(srcs, "mr_y_qc__post__aut.parquet"),
    mr_y_qc__post__aut__fmri__fov__cutoff__dorsal_max >= 65
  )

  fMRI_ventral_cutoff_score <- get_abcd_exclusion_crit(
    fs::path(srcs, "mr_y_qc__post__aut.parquet"),
    mr_y_qc__post__aut__fmri__fov__cutoff__ventral_max >= 60
  )

  Derived_results_exist <- get_abcd_exclusion_crit(
    fs::path(srcs, "mr_y_tfmri__sst__cgvfx__aseg.parquet"),
    is.na(mr_y_tfmri__sst__cgvfx__aseg__cwm__lh_beta)
  )

  purrr::reduce(
    .x = list(
      tfMRI_series_passed_rawQC,
      T1_series_passed_rawQC,
      behavior_passed,
      task_had_no_glitch,
      E_prime_timing_match_OR_ignore_E_prime_mismatch,
      fMRI_B0_unwarp_available,
      FreeSurfer_QC_not_failed,
      fMRI_manual_post_processing_QC_not_failed,
      fMRI_registration_to_T1w,
      fMRI_dorsal_cutoff_score,
      fMRI_ventral_cutoff_score,
      Derived_results_exist
    ),
    .f = dplyr::full_join,
    dplyr::join_by(sub, ses),
  ) |>
    dplyr::distinct(sub, ses) |>
    dplyr::mutate(task = "sst") |>
    dplyr::collect()
}

get_abcd_excl_nback <- function(srcs) {
  # https://docs.abcdstudy.org/latest/documentation/imaging/type_qc.html#nback-task-fmri-data-recommended-for-inclusion
  # everything except censoring

  tfMRI_series_passed_rawQC <- get_abcd_exclusion_crit(
    fs::path(srcs, "mr_y_qc__raw__tfmri__nback.parquet"),
    !(mr_y_qc__raw__tfmri__nback__pass__qc__comp_count > 0)
  )

  T1_series_passed_rawQC <- get_abcd_exclusion_crit(
    fs::path(srcs, "mr_y_qc__raw__smri__t1.parquet"),
    !(mr_y_qc__raw__smri__t1__pass__qc__comp_count > 0)
  )

  behavior_passed <- get_abcd_exclusion_crit(
    fs::path(srcs, "mr_y_tfmri__nback__beh.parquet"),
    !(mr_y_tfmri__nback__beh__qc_indicator == 1)
  )

  E_prime_timing_match_OR_ignore_E_prime_mismatch <- get_abcd_exclusion_crit(
    fs::path(srcs, "mr_y_qc__raw__tfmri__nback.parquet"),
    !(mr_y_qc__raw__tfmri__nback__eprime__match_indicator == 1 |
      mr_y_qc__raw__tfmri__nback__eprime__tdiff__ign_indicator < 1)
  )

  fMRI_B0_unwarp_available <- get_abcd_exclusion_crit(
    fs::path(srcs, "mr_y_qc__post__aut.parquet"),
    !(mr_y_qc__post__aut__fmri__b0__unwarp_indicator == 1)
  )

  FreeSurfer_QC_not_failed <- get_abcd_exclusion_crit(
    fs::path(srcs, "mr_y_qc__post__man__fsurf.parquet"),
    (mr_y_qc__post__man__fsurf_score == 0)
  )

  fMRI_manual_post_processing_QC_not_failed <- get_abcd_exclusion_crit(
    fs::path(srcs, "mr_y_qc__post__man__fmri.parquet"),
    (mr_y_qc__post__man__fmri_score == 0)
  )

  fMRI_registration_to_T1w <- get_abcd_exclusion_crit(
    fs::path(srcs, "mr_y_qc__post__aut.parquet"),
    mr_y_qc__post__aut__fmri__rigid_score >= 19
  )

  fMRI_dorsal_cutoff_score <- get_abcd_exclusion_crit(
    fs::path(srcs, "mr_y_qc__post__aut.parquet"),
    mr_y_qc__post__aut__fmri__fov__cutoff__dorsal_max >= 65
  )

  fMRI_ventral_cutoff_score <- get_abcd_exclusion_crit(
    fs::path(srcs, "mr_y_qc__post__aut.parquet"),
    mr_y_qc__post__aut__fmri__fov__cutoff__ventral_max >= 60
  )

  Derived_results_exist <- get_abcd_exclusion_crit(
    fs::path(srcs, "mr_y_tfmri__nback__0b__aseg.parquet"),
    is.na(mr_y_tfmri__nback__0b__aseg__cbwm__lh_beta)
  )

  purrr::reduce(
    .x = list(
      tfMRI_series_passed_rawQC,
      T1_series_passed_rawQC,
      behavior_passed,
      E_prime_timing_match_OR_ignore_E_prime_mismatch,
      fMRI_B0_unwarp_available,
      FreeSurfer_QC_not_failed,
      fMRI_manual_post_processing_QC_not_failed,
      fMRI_registration_to_T1w,
      fMRI_dorsal_cutoff_score,
      fMRI_ventral_cutoff_score,
      Derived_results_exist
    ),
    .f = dplyr::full_join,
    dplyr::join_by(sub, ses),
  ) |>
    dplyr::distinct(sub, ses) |>
    dplyr::mutate(task = "nback") |>
    dplyr::collect()
}

get_abcd_exclusion_official <- function(phenotypes) {
  sst <- get_abcd_excl_sst(phenotypes)
  nback <- get_abcd_excl_nback(phenotypes)
  mid <- get_abcd_excl_mid(phenotypes)
  rest <- get_abcd_excl_rest(phenotypes)
  dplyr::bind_rows(sst, nback, mid, rest)
}

get_abcd_subs_from_srcs <- function(source) {
  stringr::str_extract(source, "(?<=sub-)([[:alpha:]]|[[:digit:]])+") |>
    unique()
}

get_abcd_exclusion_demographics <- function(abcd_source, abcd_demographics) {
  entities <- get_abcd_entities(abcd_source) |> dplyr::distinct(sub, ses)

  # https://docs.abcdstudy.org/latest/documentation/release_notes/6_0.html#known-issue-scanner-issue-for-site15
  site15 <- abcd_demographics |>
    dplyr::filter(
      deviceserialnumber ==
        "d06b7afa24778649d21924836f3816031db95e5c7a5db0e475f640bb805cdf9f",
      mr_y_adm__info_dt < "2018-08-01"
    ) |>
    dplyr::distinct(sub, ses)

  dplyr::anti_join(
    entities,
    dplyr::semi_join(
      entities,
      dplyr::select(abcd_demographics, sub, ses, age, sex, bmi) |> na.omit(),
      by = dplyr::join_by(sub, ses)
    ),
    by = dplyr::join_by(sub, ses)
  ) |>
    dplyr::bind_rows(site15) |>
    dplyr::distinct()
}

exclude_bad_abcd_scan <- function(d, abcd_source, phenotypes) {
  official <- get_abcd_exclusion_official(phenotypes)
  runs <- get_abcd_entities(abcd_source) |>
    dplyr::select(-src) |>
    dplyr::filter((task == "rest" & run > 4) | (!(task == "rest") & run > 2))

  to_exclude <- dplyr::bind_rows(
    list(official = official, extra_runs = runs),
    .id = "notes"
  ) |>
    dplyr::summarise(
      notes = stringr::str_c(notes, collapse = "|"),
      .by = c(sub, ses, task, run)
    )

  d |> dplyr::anti_join(to_exclude, by = dplyr::join_by(sub, ses, task, run))
}
