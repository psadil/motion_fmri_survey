library(targets)
library(tarchetypes)
library(crew)

# TODO: convert figures to tikz

source("R/hcp.R")
source("R/ukb.R")
source("R/abcd.R")
source("R/mriqc.R")
source("R/spacetop.R")
source("R/utils.R")
source("R/spectrum.R")
source("R/qc.R")
source("R/figures.R")
source("R/cost.R")

targets::tar_option_set(
  trust_timestamps = TRUE,
  format = "parquet"
  # controller = crew::crew_controller_local(workers = 12)
)


list(
  tar_target(
    ukb_source,
    "data/motion/derivatives/ukb.parquet",
    format = "file"
  ),
  tar_target(ukb_eprime, "data/25748_2_0.txt", format = "file"),
  tar_target(ukb_responses, get_ukb_responses(ukb_eprime)),
  tar_target(hcpya_events, get_hcpya_events("data/hcp_evs", hcpya)),
  tar_target(ukb_beh, "data/ukb677207_bulk.parquet", format = "file"),
  tar_target(ukb_events, get_ukb_design()),
  tar_target(
    ukb_exclusion,
    get_ukb_exclusion(bulk = ukb_beh),
    deployment = "main"
  ),
  tar_target(ukb, get_ukb(src = ukb_source, exclusion = ukb_exclusion)),
  tar_target(
    spacetop_source,
    "data/motion/derivatives/spacetop.parquet",
    format = "file"
  ),
  tar_target(spacetop, get_spacetop(spacetop_source), deployment = "main"),
  tar_target(
    hcpya_source,
    "data/motion/derivatives/hcpya.parquet",
    format = "file"
  ),
  tar_target(
    hcpya_exclusion_src,
    "data/exclusion/hcp_ya_qc_issue_exclusion.tsv",
    format = "file"
  ),
  tar_target(
    hcpya_exclusion,
    get_hcpya_exclusion_by_other(hcpya_exclusion_src),
    deployment = "main"
  ),
  tar_target(
    hcpya,
    get_hcp(hcpya_source, hcpya_exclusion),
    deployment = "main"
  ),
  tar_target(
    hcpa_source,
    "data/motion/derivatives/hcpa.parquet",
    format = "file"
  ),
  tar_target(hcpa, get_hcp(hcpa_source, hcpa_exclusion), deployment = "main"),
  tar_target(
    hcpd_source,
    "data/motion/derivatives/hcpd.parquet",
    format = "file"
  ),
  tar_target(
    hcpd_exclusion_src,
    "data/exclusion/HCD_LS2.0_anat_anomalies_Rad_reads.xlsx",
    format = "file"
  ),
  tar_target(hcpd_exclusion, get_hcp_ls_exclusion_by_other(hcpd_exclusion_src)),
  tar_target(
    hcpa_exclusion_src,
    "data/exclusion/HCA_LS2.0_anat_anomalies_Rad_reads.xlsx",
    format = "file"
  ),
  tar_target(hcpa_exclusion, get_hcp_ls_exclusion_by_other(hcpa_exclusion_src)),
  tar_target(hcpd, get_hcp(hcpd_source, hcpd_exclusion)),
  tar_target(
    abcd_source,
    fs::dir_ls(
      "data/abcc-4-0-0/abcc-xcp_d_v0.13.0",
      glob = "*run*motion.tsv",
      recurse = TRUE
    ),
    format = "file"
  ),
  tar_target(
    abcd_phenotypes,
    "data/abcc-4-0-0/rawdata/phenotype",
    format = "file"
  ),
  tar_target(
    abcd_exclusion,
    get_abcd_exclusion_demographics(abcd_source, abcd_demographics)
  ),
  tar_target(
    abcd_events_src,
    "data/events/derivatives/abcd.parquet",
    format = "file"
  ),
  tar_target(abcd_events, get_abcd_events(abcd_events_src)),
  tar_target(abcd, get_abcd(abcd_source, abcd_exclusion, abcd_phenotypes)),
  tarchetypes::tar_group_by(
    datasets,
    bind_datasets(list(
      hcpya = hcpya,
      hcpa = hcpa,
      hcpd = hcpd,
      ukb = ukb,
      abcd = abcd,
      spacetop = spacetop
    )),
    dataset
  ),
  tar_target(
    by_time,
    summarise_by_time(datasets),
    pattern = map(datasets),
    deployment = "main"
  ),
  tar_target(
    by_run,
    summarise_by(datasets, .cols = c(dataset, task, scan, ses, sub, filtered)),
    pattern = map(datasets),
    deployment = "main"
  ),
  tar_target(ukb_demographics, get_ukb_demographics(ukb_beh)),
  tar_target(
    hcpya_unrestricted,
    "data/unrestricted_martin_2_5_2024_10_18_12.csv",
    format = "file"
  ),
  tar_target(
    hcpya_restricted,
    "data/RESTRICTED_martin_2_5_2024_10_18_28.csv",
    format = "file"
  ),
  tar_target(
    hcpya_demographics,
    get_hcpya_demographics(hcpya_unrestricted, hcpya_restricted, hcpya)
  ),
  tar_target(
    spacetop_demographics_src,
    "data/participants.tsv",
    format = "file"
  ),
  tar_target(
    spacetop_demographics,
    get_spacetop_demographics(spacetop_demographics_src)
  ),
  tar_target(hcpd_demographics, get_hcpd_demographics(hcpd)),
  tar_target(hcpa_demographics, get_hcpa_demographics(hcpa)),
  tar_target(abcd_demographics, get_abcd_demographics(abcd_source)),
  tar_target(
    demographics,
    get_demographics(
      datasets = list(
        ukb = ukb_demographics,
        hcpya = hcpya_demographics,
        spacetop = spacetop_demographics,
        hcpd = hcpd_demographics,
        hcpa = hcpa_demographics,
        abcd = abcd_demographics
      ),
      by_run = by_run
    )
  ),
  tar_target(
    demographics_tsv,
    write_demographics(demographics, "data/demographics/demographics.tsv"),
    format = "file"
  ),
  tar_target(mriqc_bold, get_mriqc_bold()),
  tar_target(
    lost_strict,
    get_lost_strict(datasets, by_run = by_run),
    pattern = map(datasets)
  ),
  tar_target(lost_lenient, get_lost_lenient(by_run = by_run)),
  tar_target(lost, bind_lost(strict = lost_strict, lenient = lost_lenient)),
  tar_target(
    hcpa_spectrum_src,
    "data/motion/derivatives/hcpa_spectrum.parquet",
    format = "file"
  ),
  tar_target(
    hcpa_spectrum,
    get_hcpa_spectrum(src = hcpa_spectrum_src, hcpa = hcpa, by_run = by_run)
  ),
  tar_target(
    hcpd_spectrum_src,
    "data/motion/derivatives/hcpd_spectrum.parquet",
    format = "file"
  ),
  tar_target(
    hcpd_spectrum,
    get_hcpd_spectrum(hcpd_spectrum_src, hcpd, by_run = by_run)
  ),
  tar_target(
    hcpya_spectrum_src,
    "data/motion/derivatives/hcpya_spectrum.parquet",
    format = "file"
  ),
  tar_target(
    hcpya_spectrum,
    get_hcpya_spectrum(src = hcpya_spectrum_src, hcpya = hcpya, by_run = by_run)
  ),
  tar_target(
    abcd_spectrum_src,
    "data/motion/derivatives/abcc_spectrum.parquet",
    format = "file"
  ),
  tar_target(
    abcd_spectrum,
    get_abcd_spectrum(src = abcd_spectrum_src, by_run = by_run)
  ),
  tar_target(
    spacetop_spectrum_src,
    "data/motion/derivatives/spacetop_spectrum.parquet",
    format = "file"
  ),
  tar_target(
    spacetop_spectrum,
    get_spacetop_spectrum(
      src = spacetop_spectrum_src,
      spacetop = spacetop,
      by_run = by_run
    )
  ),
  tar_target(
    ukb_spectrum_src,
    "data/motion/derivatives/ukb_spectrum.parquet",
    format = "file"
  ),
  tar_target(
    ukb_spectrum,
    get_ukb_spectrum(src = ukb_spectrum_src, ukb = ukb, by_run = by_run)
  ),
  tar_target(mriqc_src, "data/bold", format = "file"),
  tar_target(
    fig_age_chart,
    make_fig_demographics(by_run, demographics, mriqc_src),
    packages = c("patchwork"),
    format = "qs"
  ),
  tar_target(
    fig_by_time_hcpya,
    make_fig_by_time_hcpya(by_time, hcpya_events = hcpya_events),
    packages = c("patchwork"),
    format = "qs"
  ),
  tar_target(
    hcpd_event_times,
    "data/events_consolidated/hcp_dev_event_times.csv",
    format = "file"
  ),
  tar_target(
    fig_by_time_hcpd,
    make_fig_by_time_hcpd(
      hcpd = hcpd,
      hcpd_event_times = hcpd_event_times,
      by_time = by_time
    ),
    packages = c("patchwork"),
    format = "qs"
  ),
  tar_target(
    hcpa_event_times,
    "data/events_consolidated/hcp_aging_event_times.csv",
    format = "file"
  ),
  tar_target(
    fig_by_time_hcpa,
    make_fig_by_time_hcpa(
      hcpa = hcpa,
      hcpa_event_times = hcpa_event_times,
      by_time = by_time
    ),
    packages = c("patchwork"),
    format = "qs"
  ),
  tar_target(
    fig_by_time_ukb,
    make_fig_by_time_ukb(
      ukb_events = ukb_events,
      ukb_responses = ukb_responses,
      by_time = by_time
    ),
    packages = c("patchwork"),
    format = "qs"
  ),
  tar_target(
    abcd_events2,
    "data/events_consolidated/abcd_events.tsv",
    format = "file"
  ),
  tar_target(
    fig_by_time_abcd,
    make_fig_by_time_abcd(abcd_events = abcd_events2, by_time = by_time),
    packages = c("patchwork"),
    format = "qs"
  ),
  tar_target(
    fig_by_time_spacetop,
    make_fig_by_time_spacetop(by_time = by_time),
    packages = c("patchwork"),
    format = "qs"
  ),
  tar_target(
    fig_by_run,
    make_fig_by_run(by_run = by_run, demographics = demographics),
    packages = c("patchwork"),
    format = "qs"
  ),
  tar_target(
    fig_all_motion_exclusion,
    make_fig_all_motion_exclusion(lost = lost),
    packages = c("patchwork"),
    format = "qs"
  ),
  tar_target(
    fig_all_motion_exclusion2,
    make_fig_all_motion_exclusion(lost = lost, filtered = TRUE),
    packages = c("patchwork"),
    format = "qs"
  ),
  tar_target(group_highlow, make_group_highlow(by_run)),
  tar_target(
    fig_groups,
    make_fig_cluster(by_run, group_highlow),
    packages = c("patchwork"),
    format = "qs"
  ),
  tar_target(
    fig_lost_by_group,
    make_fig_lost_by_group(by_run, demographics),
    packages = c("patchwork"),
    format = "qs"
  ),
  tar_target(
    lost_strict2,
    get_lost_strict_for_modeling(datasets, by_run = by_run),
    pattern = map(datasets)
  ),
  tar_target(
    params,
    list(k0 = 0.3025303, k1 = 164.3939, k2 = 7128.788),
    format = "qs"
  ),
  tar_target(cost_plot, plot_cost(by_run), format = "qs"),
  tar_target(cost_lost_plot, plot_cost_lost(by_run), format = "qs"),
  tar_target(
    acc_decrease_plot,
    plot_acc_decrease(by_run, params = params),
    format = "qs"
  ),
  tar_target(
    fig_cost,
    make_fig_cost(cost_plot, cost_lost_plot, acc_decrease_plot),
    format = "qs",
    packages = "patchwork"
  ),
  tar_target(key_measures, get_key_measures(), format = "qs"),
  tar_target(performance_src, "data/performance", format = "file"),
  tar_target(performance, get_performance(performance_src)),
  tar_target(
    fig_performance,
    make_fig_performance(performance, key_measures),
    packages = "patchwork",
    format = "qs"
  ),
  tar_target(fig_power, make_fig_power(by_run, lost), format = "qs"),
  tar_target(fig_hcp_psych, make_fig_hcp_psych(by_run), format = "qs"),
  tar_target(timeseries_src, "data/timeseries", format = "file"),
  tar_target(
    ukb_subs,
    get_ukb_subs(by_run, timeseries_src, n = 500),
    format = "qs"
  ),
  tar_target(
    qc_fd_ukb_real0,
    get_qc_fd_ukb_real(ukb, ukb_subs),
    deployment = "main"
  ),
  tarchetypes::tar_group_by(qc_fd_ukb_real, qc_fd_ukb_real0, sub),
  tar_target(iter, seq_len(10), format = "qs"),
  tar_target(
    qc_fd_ukb0,
    get_qc_fd_ukb(dplyr::select(qc_fd_ukb_real, -tar_group), new_iter = iter),
    pattern = cross(qc_fd_ukb_real, iter)
  ),
  tar_target(type_id, c("clean", "raw", "rawraw"), format = "qs"), # for ukb timeseries
  tarchetypes::tar_group_by(
    qc_fd_ukb,
    tidyr::crossing(qc_fd_ukb0, threshold = seq(0.01, .2, by = 0.02)),
    threshold,
    sub
  ),
  tar_target(
    mac_prep,
    get_mac_prep(
      d = dplyr::select(qc_fd_ukb, -tar_group),
      timeseries_src = timeseries_src,
      type_id = type_id
    ),
    pattern = cross(qc_fd_ukb, type_id)
  ),
  tar_target(
    mac0,
    get_mac(mac_prep, real = qc_fd_ukb_real, timeseries_src = timeseries_src),
    pattern = map(mac_prep)
  ),
  tar_target(
    mac,
    mac0 |>
      dplyr::summarise(
        v = var(mac),
        sem = sd(mac) / sqrt(dplyr::n()),
        mac = mean(mac),
        .by = c(threshold, type, filtered)
      ),
  ),
  tar_target(fig_mac, make_fig_mac(mac), format = "qs"),
  tar_target(
    orig_counts,
    get_orig_counts(
      ukb_source = ukb_source,
      spacetop_source = spacetop_source,
      abcd_source = abcd_source,
      hcpya_source = hcpya_source,
      hcpa_source = hcpa_source,
      hcpd_source = hcpd_source
    )
  ),
  tar_target(
    freq_limits,
    prep_tbl_freq_limits(
      hcpa_spectrum = hcpa_spectrum,
      hcpd_spectrum = hcpd_spectrum,
      hcpya_spectrum = hcpya_spectrum,
      spacetop_spectrum = spacetop_spectrum,
      ukb_spectrum = ukb_spectrum,
      abcd_spectrum = abcd_spectrum,
      demographics = demographics
    )
  ),
  tar_target(fig_hcpa_spectrum, plot_spectrum(hcpa_spectrum), format = "qs"),
  tar_target(fig_hcpd_spectrum, plot_spectrum(hcpd_spectrum), format = "qs"),
  tar_target(fig_hcpya_spectrum, plot_spectrum(hcpya_spectrum), format = "qs"),
  tar_target(fig_ukb_spectrum, plot_spectrum(ukb_spectrum), format = "qs"),
  tar_target(
    fig_spacetop_spectrum,
    plot_spectrum(spacetop_spectrum),
    format = "qs"
  ),
  tar_target(
    fig_abcd_spectrum_00,
    plot_spectrum(dplyr::filter(abcd_spectrum, ses == "Baseline")),
    format = "qs"
  ),
  tar_target(
    fig_abcd_spectrum_02,
    plot_spectrum(dplyr::filter(abcd_spectrum, ses == "Year2")),
    format = "qs"
  ),
  tar_target(
    fig_abcd_spectrum_04,
    plot_spectrum(dplyr::filter(abcd_spectrum, ses == "Year4")),
    format = "qs"
  ),
  tar_target(
    fig_abcd_spectrum_06,
    plot_spectrum(dplyr::filter(abcd_spectrum, ses == "Year6")),
    format = "qs"
  ),
  tar_target(icc, get_icc(by_run)),
  tar_target(
    bpm_src,
    "data/abcc-4-0-0/rawdata/phenotype/mh_t_bpm.parquet",
    format = "file"
  ),
  tar_target(
    compare_datasets_fit,
    get_compare_datasets_fit(by_run, demographics, bpm_src),
    format = "qs"
  ),
  tar_target(rest_fit, get_rest_fit(by_run, demographics), format = "qs"),
  tar_target(
    mr_y_qc__mot,
    "data/abcc-4-0-0/rawdata/phenotype/mr_y_qc__mot.parquet",
    format = "file"
  ),
  tar_target(fig_abcc_abcd, plot_abcc_abcd(abcd, mr_y_qc__mot), format = "qs"),
  tarchetypes::tar_quarto(report, quiet = FALSE)
)
