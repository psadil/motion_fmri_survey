get_abcd <- function(sources, abcd_exclusion) {
  to_exclude <- abcd_exclusion |>
    arrow::as_arrow_table() |>
    dplyr::mutate(
      sub = arrow::cast(sub, arrow::large_utf8()),
      ses = arrow::cast(ses, arrow::large_utf8()),
      task = arrow::cast(task, arrow::large_utf8()),
      run = arrow::cast(run, arrow::large_utf8())
    )

  arrow::open_dataset(sources) |>
    dplyr::anti_join(to_exclude) |>
    dplyr::filter(t > 0) |>
    dplyr::mutate(
      time = t * 0.8,
      run = as.integer(run),
      ped = "AP",
      scan = as.integer(run)
    ) |>
    dplyr::collect() |>
    add_ped() |>
    do_casting()
}

get_abcd_events <- function(abcd_design = "/Users/psadil/data/events/derivatives/abcd.parquet") {
  events <- arrow::open_dataset(abcd_design)

  cue_design <- events |>
    dplyr::filter(task == "nback") |>
    dplyr::filter(trial_type == "cue") |>
    dplyr::arrange(onset) |>
    dplyr::collect() |>
    dplyr::mutate(
      event = 1:dplyr::n(),
      .by = c(run, task, sub, ses)
    ) |>
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
    dplyr::mutate(
      event = 1:dplyr::n(),
      .by = c(run, task, sub, ses)
    ) |>
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
    dplyr::mutate(
      event = 1:dplyr::n(),
      .by = c(run, task, sub, ses)
    ) |>
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
    dplyr::mutate(
      event = 1:dplyr::n(),
      .by = c(run, task, sub, ses)
    ) |>
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


get_abcds_with_scanner <- function(.data) {
  with_na <- .data |>
    dplyr::group_nest(sub, ses) |>
    dplyr::mutate(
      anyna = purrr::map_lgl(
        data,
        ~ any(is.na(.x$deviceserialnumber)) && !all(is.na(.x$deviceserialnumber))
      )
    ) |>
    dplyr::filter(anyna) |>
    tidyr::unnest(data) |>
    na.omit() |>
    dplyr::select(-anyna)

  .data |>
    dplyr::anti_join(with_na) |> # de-duplicate
    dplyr::bind_rows(with_na) # add the correct entries back
}


get_abcd_demographics <- function(source) {
  abcd_subs <- arrow::open_dataset(source) |>
    dplyr::distinct(sub, ses) |>
    dplyr::collect()

  main <- read_nda("data/demographics/image03.txt") |>
    dplyr::rename(ses = visit, age = interview_age) |>
    dplyr::mutate(
      sub = stringr::str_remove(src_subject_id, "_"),
      interview_date = lubridate::as_date(interview_date, format = "%m/%d/%Y")
    ) |>
    dplyr::semi_join(dplyr::select(abcd_subs, sub)) |>
    dplyr::mutate(
      ses = dplyr::case_match(
        ses,
        "baseline_year_1_arm_1" ~ "baselineYear1Arm1",
        "2_year_follow_up_y_arm_1" ~ "2YearFollowUpYArm1",
        "4_year_follow_up_y_arm_1" ~ "4YearFollowUpYArm1",
        .default = ses
      ),
      age = age / 12,
      sex = dplyr::case_match(
        sex,
        "F" ~ "Female",
        "M" ~ "Male",
        .default = sex
      )
    ) |>
    dplyr::distinct(sub, deviceserialnumber, interview_date, age, sex, ses) |>
    get_abcds_with_scanner() |>
    dplyr::slice_min(order_by = interview_date, n = 1, by = c(sub, ses)) |>
    dplyr::slice_min(order_by = age, n = 1, by = c(sub, ses)) |> # errors in dataset
    dplyr::select(-interview_date)

  age <- readr::read_csv("data/abcd_y_lt.csv") |>
    dplyr::rename(sub = src_subject_id, ses = eventname) |>
    dplyr::mutate(
      sub = stringr::str_remove(sub, "_"),
      ses = dplyr::case_match(
        ses,
        "baseline_year_1_arm_1" ~ "baselineYear1Arm1",
        "2_year_follow_up_y_arm_1" ~ "2YearFollowUpYArm1",
        "4_year_follow_up_y_arm_1" ~ "4YearFollowUpYArm1",
        .default = ses
      ),
      age = interview_age / 12
    ) |>
    dplyr::filter(!is.na(age)) |>
    dplyr::select(sub, ses, age)

  sex <- readr::read_csv("data/gish_y_gi.csv") |>
    dplyr::rename(sub = src_subject_id, ses = eventname) |>
    dplyr::mutate(
      sub = stringr::str_remove(sub, "_"),
      ses = dplyr::case_match(
        ses,
        "baseline_year_1_arm_1" ~ "baselineYear1Arm1",
        "2_year_follow_up_y_arm_1" ~ "2YearFollowUpYArm1",
        "4_year_follow_up_y_arm_1" ~ "4YearFollowUpYArm1",
        .default = ses
      ),
      sex = kbi_sex_assigned_at_birth
    ) |>
    dplyr::select(sub, ses, sex) |>
    tidyr::pivot_wider(names_from = ses, values_from = sex) |>
    dplyr::mutate(
      sex = dplyr::if_else(is.na(baselineYear1Arm1), `1_year_follow_up_y_arm_1`, baselineYear1Arm1),
      sex = dplyr::if_else(is.na(sex), `2YearFollowUpYArm1`, sex),
      sex = dplyr::if_else(is.na(sex), `3_year_follow_up_y_arm_1`, sex),
      sex = dplyr::if_else(is.na(sex), `4YearFollowUpYArm1`, sex),
      sex = dplyr::case_match(
        sex,
        1 ~ "Male",
        2 ~ "Female",
        .default = NA
      )
    ) |>
    dplyr::filter(!is.na(sex)) |>
    dplyr::select(sub, sex)


  bmi <- readr::read_csv("data/ph_y_anthro.csv") |>
    dplyr::rename(sub = src_subject_id, ses = eventname) |>
    dplyr::mutate(
      sub = stringr::str_remove(sub, "_"),
      ses = dplyr::case_match(
        ses,
        "baseline_year_1_arm_1" ~ "baselineYear1Arm1",
        "2_year_follow_up_y_arm_1" ~ "2YearFollowUpYArm1",
        "4_year_follow_up_y_arm_1" ~ "4YearFollowUpYArm1",
        .default = ses
      )
    ) |>
    dplyr::filter(!is.na(anthroheightcalc), !is.na(anthroweightcalc)) |>
    dplyr::mutate(
      weight = anthroweightcalc * 0.4535924,
      height = anthroheightcalc * 0.0254,
      bmi = weight / height^2
    ) |>
    dplyr::select(sub, ses, bmi) |>
    dplyr::filter(!is.na(bmi), bmi < 500)

  dplyr::full_join(main, bmi, by = dplyr::join_by(sub, ses)) |>
    dplyr::full_join(sex, by = dplyr::join_by(sub)) |>
    dplyr::full_join(age, by = dplyr::join_by(sub, ses)) |>
    dplyr::filter(ses %in% c("baselineYear1Arm1", "2YearFollowUpYArm1", "4YearFollowUpYArm1")) |>
    dplyr::mutate(
      age = dplyr::if_else(is.na(age.y), age.x, age.y),
      sex = dplyr::if_else(is.na(sex.y), sex.x, sex.y),
    ) |>
    dplyr::select(sub, ses, age, sex, deviceserialnumber, bmi)
}


get_abcd_excl_rest <- function() {
  # https://wiki.abcdstudy.org/release-notes/imaging/quality-control.html#rs-fmri-data-recommended-for-inclusion
  # everything except censoring
  srcs <- here::here("data", "abcd_tabular", "imaging")
  rsfMRI_tfMRI_series_passed_rawQC <- readr::read_csv(
    fs::path(srcs, "mri_y_qc_raw_rsfmr.csv"),
    col_select = c("sub" = "src_subject_id", "ses" = "eventname", "iqc_rsfmri_ok_ser")
  ) |>
    dplyr::filter(!(iqc_rsfmri_ok_ser > 0))

  T1_series_passed_rawQC <- readr::read_csv(
    fs::path(srcs, "mri_y_qc_raw_smr_t1.csv"),
    col_select = c("sub" = "src_subject_id", "ses" = "eventname", "iqc_t1_ok_ser")
  ) |>
    dplyr::filter(!(iqc_t1_ok_ser > 0))

  fMRI_B0_unwarp_available <- readr::read_csv(
    fs::path(srcs, "mri_y_qc_auto_post.csv"),
    col_select = c("sub" = "src_subject_id", "ses" = "eventname", "apqc_fmri_bounwarp_flag")
  ) |>
    dplyr::filter(!(apqc_fmri_bounwarp_flag == 1))

  FreeSurfer_QC_not_failed <- readr::read_csv(
    fs::path(srcs, "mri_y_qc_man_fsurf.csv"),
    col_select = c("sub" = "src_subject_id", "ses" = "eventname", "fsqc_qc")
  ) |>
    dplyr::filter(fsqc_qc == 0)

  fMRI_manual_post_processing_QC_not_failed <- readr::read_csv(
    fs::path(srcs, "mri_y_qc_man_post_fmr.csv"),
    col_select = c("sub" = "src_subject_id", "ses" = "eventname", "fmri_postqc_qc")
  ) |>
    dplyr::filter(fmri_postqc_qc == 0)

  fMRI_registration_to_T1w <- readr::read_csv(
    fs::path(srcs, "mri_y_qc_auto_post.csv"),
    col_select = c("sub" = "src_subject_id", "ses" = "eventname", "apqc_fmri_regt1_rigid")
  ) |>
    dplyr::filter(!(apqc_fmri_regt1_rigid < 19))

  fMRI_dorsal_cutoff_score <- readr::read_csv(
    fs::path(srcs, "mri_y_qc_auto_post.csv"),
    col_select = c("sub" = "src_subject_id", "ses" = "eventname", "apqc_fmri_fov_cutoff_dorsal")
  ) |>
    dplyr::filter(!(apqc_fmri_fov_cutoff_dorsal < 65))

  fMRI_ventral_cutoff_score <- readr::read_csv(
    fs::path(srcs, "mri_y_qc_auto_post.csv"),
    col_select = c("sub" = "src_subject_id", "ses" = "eventname", "apqc_fmri_fov_cutoff_ventral")
  ) |>
    dplyr::filter(!(apqc_fmri_fov_cutoff_ventral < 60))

  Derived_results_exist <- readr::read_csv(
    fs::path(srcs, "mri_y_rsfmr_cor_gp_gp.csv"),
    col_select = c("sub" = "src_subject_id", "ses" = "eventname", "rsfmri_c_ngd_dt_ngd_sa")
  ) |>
    dplyr::filter(is.na(rsfmri_c_ngd_dt_ngd_sa))

  purrr::reduce(
    .x = list(
      rsfMRI_tfMRI_series_passed_rawQC,
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
    dplyr::mutate(task = "rest")
}


get_abcd_excl_mid <- function() {
  # https://wiki.abcdstudy.org/release-notes/imaging/quality-control.html#rs-fmri-data-recommended-for-inclusion
  # everything except censoring
  srcs <- here::here("data", "abcd_tabular", "imaging")
  MID_tfMRI_series_passed_rawQC <- readr::read_csv(
    fs::path(srcs, "mri_y_qc_raw_tfmr_mid.csv"),
    col_select = c("sub" = "src_subject_id", "ses" = "eventname", "iqc_mid_ok_ser")
  ) |>
    dplyr::filter(!(iqc_mid_ok_ser > 0))

  T1_series_passed_rawQC <- readr::read_csv(
    fs::path(srcs, "mri_y_qc_raw_smr_t1.csv"),
    col_select = c("sub" = "src_subject_id", "ses" = "eventname", "iqc_t1_ok_ser")
  ) |>
    dplyr::filter(!(iqc_t1_ok_ser > 0))

  MID_behavior_passed <- readr::read_csv(
    fs::path(srcs, "mri_y_tfmr_mid_beh.csv"),
    col_select = c("sub" = "src_subject_id", "ses" = "eventname", "tfmri_mid_beh_performflag")
  ) |>
    dplyr::filter(!(tfmri_mid_beh_performflag == 1))

  MID_E_prime_timing_match_OR_ignore_E_prime_mismatch <- readr::read_csv(
    fs::path(srcs, "mri_y_qc_raw_tfmr_mid.csv"),
    col_select = c("sub" = "src_subject_id", "ses" = "eventname", "iqc_mid_ep_t_series_match", "eprime_mismatch_ok_mid")
  ) |>
    dplyr::filter(!(iqc_mid_ep_t_series_match == 1 | eprime_mismatch_ok_mid == 1))

  fMRI_B0_unwarp_available <- readr::read_csv(
    fs::path(srcs, "mri_y_qc_auto_post.csv"),
    col_select = c("sub" = "src_subject_id", "ses" = "eventname", "apqc_fmri_bounwarp_flag")
  ) |>
    dplyr::filter(!(apqc_fmri_bounwarp_flag == 1))

  FreeSurfer_QC_not_failed <- readr::read_csv(
    fs::path(srcs, "mri_y_qc_man_fsurf.csv"),
    col_select = c("sub" = "src_subject_id", "ses" = "eventname", "fsqc_qc")
  ) |>
    dplyr::filter(fsqc_qc == 0)

  fMRI_manual_post_processing_QC_not_failed <- readr::read_csv(
    fs::path(srcs, "mri_y_qc_man_post_fmr.csv"),
    col_select = c("sub" = "src_subject_id", "ses" = "eventname", "fmri_postqc_qc")
  ) |>
    dplyr::filter(fmri_postqc_qc == 0)

  fMRI_registration_to_T1w <- readr::read_csv(
    fs::path(srcs, "mri_y_qc_auto_post.csv"),
    col_select = c("sub" = "src_subject_id", "ses" = "eventname", "apqc_fmri_regt1_rigid")
  ) |>
    dplyr::filter(!(apqc_fmri_regt1_rigid < 19))

  fMRI_dorsal_cutoff_score <- readr::read_csv(
    fs::path(srcs, "mri_y_qc_auto_post.csv"),
    col_select = c("sub" = "src_subject_id", "ses" = "eventname", "apqc_fmri_fov_cutoff_dorsal")
  ) |>
    dplyr::filter(!(apqc_fmri_fov_cutoff_dorsal < 65))

  fMRI_ventral_cutoff_score <- readr::read_csv(
    fs::path(srcs, "mri_y_qc_auto_post.csv"),
    col_select = c("sub" = "src_subject_id", "ses" = "eventname", "apqc_fmri_fov_cutoff_ventral")
  ) |>
    dplyr::filter(!(apqc_fmri_fov_cutoff_ventral < 60))

  Derived_results_exist <- readr::read_csv(
    fs::path(srcs, "mri_y_tfmr_mid_arvn_aseg.csv"),
    col_select = c("sub" = "src_subject_id", "ses" = "eventname", "tfmri_ma_acdn_b_scs_cbwmlh")
  ) |>
    dplyr::filter(is.na(tfmri_ma_acdn_b_scs_cbwmlh))

  purrr::reduce(
    .x = list(
      MID_tfMRI_series_passed_rawQC,
      T1_series_passed_rawQC,
      MID_behavior_passed,
      MID_E_prime_timing_match_OR_ignore_E_prime_mismatch,
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
    dplyr::mutate(task = "mid")
}

get_abcd_excl_sst <- function() {
  # https://wiki.abcdstudy.org/release-notes/imaging/quality-control.html#rs-fmri-data-recommended-for-inclusion
  # everything except censoring
  srcs <- here::here("data", "abcd_tabular", "imaging")
  SST_tfMRI_series_passed_rawQC <- readr::read_csv(
    fs::path(srcs, "mri_y_qc_raw_tfmr_sst.csv"),
    col_select = c("sub" = "src_subject_id", "ses" = "eventname", "iqc_sst_ok_ser")
  ) |>
    dplyr::filter(!(iqc_sst_ok_ser > 0))

  T1_series_passed_rawQC <- readr::read_csv(
    fs::path(srcs, "mri_y_qc_raw_smr_t1.csv"),
    col_select = c("sub" = "src_subject_id", "ses" = "eventname", "iqc_t1_ok_ser")
  ) |>
    dplyr::filter(!(iqc_t1_ok_ser > 0))

  SST_behavior_passed <- readr::read_csv(
    fs::path(srcs, "mri_y_tfmr_mid_beh.csv"),
    col_select = c("sub" = "src_subject_id", "ses" = "eventname", "tfmri_mid_beh_performflag")
  ) |>
    dplyr::filter(!(tfmri_mid_beh_performflag == 1))

  SST_task_had_no_glitch <- readr::read_csv(
    fs::path(srcs, "mri_y_tfmr_sst_beh.csv"),
    col_select = c("sub" = "src_subject_id", "ses" = "eventname", "tfmri_sst_beh_glitchflag")
  ) |>
    dplyr::filter(!(tfmri_sst_beh_glitchflag == 0))

  fMRI_B0_unwarp_available <- readr::read_csv(
    fs::path(srcs, "mri_y_qc_auto_post.csv"),
    col_select = c("sub" = "src_subject_id", "ses" = "eventname", "apqc_fmri_bounwarp_flag")
  ) |>
    dplyr::filter(!(apqc_fmri_bounwarp_flag == 1))

  SST_E_prime_timing_match_OR_ignore_E_prime_mismatch <- readr::read_csv(
    fs::path(srcs, "mri_y_qc_raw_tfmr_sst.csv"),
    col_select = c("sub" = "src_subject_id", "ses" = "eventname", "iqc_sst_ep_t_series_match", "eprime_mismatch_ok_sst")
  ) |>
    dplyr::filter(!(iqc_sst_ep_t_series_match == 1 | eprime_mismatch_ok_sst == 1))

  FreeSurfer_QC_not_failed <- readr::read_csv(
    fs::path(srcs, "mri_y_qc_man_fsurf.csv"),
    col_select = c("sub" = "src_subject_id", "ses" = "eventname", "fsqc_qc")
  ) |>
    dplyr::filter(fsqc_qc == 0)

  fMRI_manual_post_processing_QC_not_failed <- readr::read_csv(
    fs::path(srcs, "mri_y_qc_man_post_fmr.csv"),
    col_select = c("sub" = "src_subject_id", "ses" = "eventname", "fmri_postqc_qc")
  ) |>
    dplyr::filter(fmri_postqc_qc == 0)

  fMRI_registration_to_T1w <- readr::read_csv(
    fs::path(srcs, "mri_y_qc_auto_post.csv"),
    col_select = c("sub" = "src_subject_id", "ses" = "eventname", "apqc_fmri_regt1_rigid")
  ) |>
    dplyr::filter(!(apqc_fmri_regt1_rigid < 19))

  fMRI_dorsal_cutoff_score <- readr::read_csv(
    fs::path(srcs, "mri_y_qc_auto_post.csv"),
    col_select = c("sub" = "src_subject_id", "ses" = "eventname", "apqc_fmri_fov_cutoff_dorsal")
  ) |>
    dplyr::filter(!(apqc_fmri_fov_cutoff_dorsal < 65))

  fMRI_ventral_cutoff_score <- readr::read_csv(
    fs::path(srcs, "mri_y_qc_auto_post.csv"),
    col_select = c("sub" = "src_subject_id", "ses" = "eventname", "apqc_fmri_fov_cutoff_ventral")
  ) |>
    dplyr::filter(!(apqc_fmri_fov_cutoff_ventral < 60))

  Derived_results_exist <- readr::read_csv(
    fs::path(srcs, "mri_y_tfmr_sst_cgvfx_aseg.csv"),
    col_select = c("sub" = "src_subject_id", "ses" = "eventname", "tfmri_sacgvf_bscs_cbwmlh")
  ) |>
    dplyr::filter(is.na(tfmri_sacgvf_bscs_cbwmlh))

  purrr::reduce(
    .x = list(
      SST_tfMRI_series_passed_rawQC,
      T1_series_passed_rawQC,
      SST_behavior_passed,
      SST_task_had_no_glitch,
      SST_E_prime_timing_match_OR_ignore_E_prime_mismatch,
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
    dplyr::mutate(task = "sst")
}

get_abcd_excl_nback <- function() {
  # https://wiki.abcdstudy.org/release-notes/imaging/quality-control.html#rs-fmri-data-recommended-for-inclusion
  # everything except censoring
  srcs <- here::here("data", "abcd_tabular", "imaging")
  nBack_tfMRI_series_passed_rawQC <- readr::read_csv(
    fs::path(srcs, "mri_y_qc_raw_tfmr_nback.csv"),
    col_select = c("sub" = "src_subject_id", "ses" = "eventname", "iqc_nback_ok_ser")
  ) |>
    dplyr::filter(!(iqc_nback_ok_ser > 0))

  T1_series_passed_rawQC <- readr::read_csv(
    fs::path(srcs, "mri_y_qc_raw_smr_t1.csv"),
    col_select = c("sub" = "src_subject_id", "ses" = "eventname", "iqc_t1_ok_ser")
  ) |>
    dplyr::filter(!(iqc_t1_ok_ser > 0))

  nBack_behavior_passed <- readr::read_csv(
    fs::path(srcs, "mri_y_tfmr_nback_beh.csv"),
    col_select = c("sub" = "src_subject_id", "ses" = "eventname", "tfmri_nback_beh_performflag")
  ) |>
    dplyr::filter(!(tfmri_nback_beh_performflag == 1))

  nBack_E_prime_timing_match_OR_ignore_E_prime_mismatch <- readr::read_csv(
    fs::path(srcs, "mri_y_qc_raw_tfmr_nback.csv"),
    col_select = c("sub" = "src_subject_id", "ses" = "eventname", "iqc_nback_ep_t_series_match", "eprime_mismatch_ok_nback")
  ) |>
    dplyr::filter(!(iqc_nback_ep_t_series_match == 1 | eprime_mismatch_ok_nback == 1))

  fMRI_B0_unwarp_available <- readr::read_csv(
    fs::path(srcs, "mri_y_qc_auto_post.csv"),
    col_select = c("sub" = "src_subject_id", "ses" = "eventname", "apqc_fmri_bounwarp_flag")
  ) |>
    dplyr::filter(!(apqc_fmri_bounwarp_flag == 1))

  FreeSurfer_QC_not_failed <- readr::read_csv(
    fs::path(srcs, "mri_y_qc_man_fsurf.csv"),
    col_select = c("sub" = "src_subject_id", "ses" = "eventname", "fsqc_qc")
  ) |>
    dplyr::filter(fsqc_qc == 0)

  fMRI_manual_post_processing_QC_not_failed <- readr::read_csv(
    fs::path(srcs, "mri_y_qc_man_post_fmr.csv"),
    col_select = c("sub" = "src_subject_id", "ses" = "eventname", "fmri_postqc_qc")
  ) |>
    dplyr::filter(fmri_postqc_qc == 0)

  fMRI_registration_to_T1w <- readr::read_csv(
    fs::path(srcs, "mri_y_qc_auto_post.csv"),
    col_select = c("sub" = "src_subject_id", "ses" = "eventname", "apqc_fmri_regt1_rigid")
  ) |>
    dplyr::filter(!(apqc_fmri_regt1_rigid < 19))

  fMRI_dorsal_cutoff_score <- readr::read_csv(
    fs::path(srcs, "mri_y_qc_auto_post.csv"),
    col_select = c("sub" = "src_subject_id", "ses" = "eventname", "apqc_fmri_fov_cutoff_dorsal")
  ) |>
    dplyr::filter(!(apqc_fmri_fov_cutoff_dorsal < 65))

  fMRI_ventral_cutoff_score <- readr::read_csv(
    fs::path(srcs, "mri_y_qc_auto_post.csv"),
    col_select = c("sub" = "src_subject_id", "ses" = "eventname", "apqc_fmri_fov_cutoff_ventral")
  ) |>
    dplyr::filter(!(apqc_fmri_fov_cutoff_ventral < 60))

  Derived_results_exist <- readr::read_csv(
    fs::path(srcs, "mri_y_tfmr_nback_0b_aseg.csv"),
    col_select = c("sub" = "src_subject_id", "ses" = "eventname", "tfmri_nback_all_4")
  ) |>
    dplyr::filter(is.na(tfmri_nback_all_4))

  purrr::reduce(
    .x = list(
      nBack_tfMRI_series_passed_rawQC,
      T1_series_passed_rawQC,
      nBack_behavior_passed,
      nBack_E_prime_timing_match_OR_ignore_E_prime_mismatch,
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
    dplyr::mutate(task = "nback")
}

get_abcd_too_short <- function(sources) {
  # can't go by max n_tr. Need to go by most common
  expected <- arrow::open_dataset(sources) |>
    dplyr::summarise(n_tr = max(t), .by = c(sub, task, ses, run)) |>
    dplyr::count(n_tr, task, ses, run) |>
    dplyr::collect() |>
    dplyr::slice_max(order_by = n, n = 1, with_ties = FALSE, by = c(task, ses, run)) |>
    dplyr::select(-n)

  arrow::open_dataset(sources) |>
    dplyr::summarise(n_tr = max(t), .by = c(sub, task, ses, run)) |>
    dplyr::collect() |>
    dplyr::anti_join(expected) |>
    dplyr::select(sub, task, ses, run) |>
    dplyr::collect()
}

get_abcd_exclusion_official <- function() {
  sst <- get_abcd_excl_sst()
  nback <- get_abcd_excl_nback()
  mid <- get_abcd_excl_mid()
  rest <- get_abcd_excl_rest()
  dplyr::bind_rows(sst, nback, mid, rest) |>
    dplyr::mutate(
      ses = dplyr::case_match(
        ses,
        "baseline_year_1_arm_1" ~ "baselineYear1Arm1",
        "2_year_follow_up_y_arm_1" ~ "2YearFollowUpYArm1",
        "4_year_follow_up_y_arm_1" ~ "4YearFollowUpYArm1",
      ),
      sub = stringr::str_remove(sub, "_")
    )
}

get_abcd_exclusion_run <- function(sources) {
  arrow::open_dataset(sources) |>
    dplyr::distinct(sub, task, ses, run) |>
    dplyr::filter(run %in% c("05", "06")) |>
    dplyr::collect()
}

get_abcd_exclusion_demographics <- function(source, demographics) {
  abcd_subs <- arrow::open_dataset(source) |>
    dplyr::distinct(sub, ses, task, run) |>
    dplyr::collect()

  dplyr::anti_join(
    abcd_subs,
    dplyr::semi_join(abcd_subs, na.omit(dplyr::select(demographics, sub, ses, age, sex, bmi)), by = dplyr::join_by(sub, ses)),
    by = dplyr::join_by(sub, ses)
  )
}


get_abcd_exclusion <- function(source, demographics) {
  all_runs <- arrow::open_dataset(source) |>
    dplyr::distinct(sub, task, ses, run) |>
    dplyr::collect()
  too_short <- all_runs |>
    dplyr::semi_join(
      get_abcd_too_short(source),
      by = dplyr::join_by(sub, task, ses, run)
    )
  official <- all_runs |>
    dplyr::semi_join(
      get_abcd_exclusion_official(),
      by = dplyr::join_by(sub, task, ses)
    )
  runs <- all_runs |>
    dplyr::semi_join(
      get_abcd_exclusion_run(source),
      by = dplyr::join_by(sub, task, ses, run)
    )

  from_demographics <- get_abcd_exclusion_demographics(source, demographics)
  dplyr::bind_rows(
    list(
      atypical_length = too_short,
      official = official,
      extra_runs = runs,
      missing_demographics = from_demographics
    ),
    .id = "notes"
  ) |>
    dplyr::summarise(
      notes = stringr::str_c(notes, collapse = "|"),
      .by = c(sub, ses, task, run)
    )

  # dplyr::bind_rows(too_short, runs) |>
  #   dplyr::distinct()
}
