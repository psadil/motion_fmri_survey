library(targets)
library(crew)

source("R/hcp.R")
source("R/ukb.R")
source("R/abcd.R")
source("R/mriqc.R")
source("R/spacetop.R")
source("R/utils.R")
source("R/spectrum.R")
source("R/qc.R")
source("R/figures.R")

targets::tar_option_set(
  trust_timestamps = TRUE,
  format = "parquet"
)


list(
  tar_target(
    ukb_source,
    "data/motion/derivatives/ukb.parquet",
    format = "file"
  ),
  tar_target(
    ukb_exclusion_src,
    "data/exclusion/ukb_exclusion.tsv",
    format = "file"
  ),
  tar_target(
    ukb_exclusion,
    get_ukb_exclusion(ukb_exclusion_src),
    deployment = "main"
  ),
  tar_target(ukb_eprime, "data/25748_2_0.txt", format = "file"),
  tar_target(ukb_responses, get_ukb_responses(ukb_eprime)),
  tar_target(hcpya_events, get_hcpya_events("data/hcp_evs", hcpya)),
  tar_target(ukb, get_ukb(ukb_source, ukb_exclusion)),
  tar_target(ukb_beh, "data/ukb677207_bulk.parquet", format = "file"),
  tar_target(ukb_events, get_ukb_design()),
  tar_target(
    spacetop_source,
    "data/motion/derivatives/spacetop.parquet",
    format = "file"
  ),
  tar_target(
    spacetop_exclusion_src,
    "data/exclusion/Spacetop_exclusions.rds",
    format = "file"
  ),
  tar_target(
    spacetop_exclusion,
    get_spacetop_exclusion(spacetop_exclusion_src),
    deployment = "main"
  ),
  tar_target(
    spacetop,
    get_spacetop(spacetop_source, spacetop_exclusion),
    deployment = "main"
  ),
  tar_target(
    hcpya_source,
    "data/motion/derivatives/hcpya.parquet",
    format = "file"
  ),
  tar_target(hcpya_exclusion, get_hcpya_exclusion(hcpya_source)),
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
  tar_target(hcpa, get_hcp(hcpa_source), deployment = "main"),
  tar_target(
    hcpd_source,
    "data/motion/derivatives/hcpd.parquet",
    format = "file"
  ),
  tar_target(
    hcpd_exclusion_src,
    "data/exclusion/hcp_dev.csv",
    format = "file"
  ),
  tar_target(hcpd_exclusion, get_hcpd_exclusion(hcpd_exclusion_src)),
  tar_target(hcpd, get_hcp(hcpd_source, hcpd_exclusion)),
  tar_target(
    abcd_source,
    "data/motion/derivatives/abcd.parquet",
    format = "file"
  ),
  tar_target(
    abcd_exclusion,
    get_abcd_exclusion(abcd_source, abcd_demographics)
  ),
  tar_target(
    abcd_events_src,
    "data/events/derivatives/abcd.parquet",
    format = "file"
  ),
  tar_target(abcd_events, get_abcd_events(abcd_events_src)),
  tar_target(abcd, get_abcd(abcd_source, abcd_exclusion)),
  tarchetypes::tar_group_by(
    datasets,
    bind_datasets(
      list(
        hcpya = hcpya,
        hcpa = hcpa,
        hcpd = hcpd,
        ukb = ukb,
        abcd = abcd,
        spacetop = spacetop
      )
    ),
    dataset,
    deployment = "main"
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
  tar_target(ukb_demographics, get_ukb_demographics(ukb_beh, ukb_source)),
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
    get_hcpya_demographics(hcpya_unrestricted, hcpya_restricted)
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
  tar_target(hcpdev_demographics, get_hcp_dev_demographics()),
  tar_target(hcpaging_demographics, get_hcp_aging_demographics()),
  tar_target(abcd_demographics, get_abcd_demographics(abcd_source)),
  tar_target(
    demographics,
    get_demographics(
      ukb = ukb_demographics,
      hcpya = hcpya_demographics,
      spacetop = spacetop_demographics,
      hcpdev = hcpdev_demographics,
      hcpaging = hcpaging_demographics,
      abcd = abcd_demographics
    ),
    format = "parquet"
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
    get_hcpa_spectrum(src = hcpa_spectrum_src, by_run = by_run)
  ),
  tar_target(
    hcpd_spectrum_src,
    "data/motion/derivatives/hcpd_spectrum.parquet",
    format = "file"
  ),
  tar_target(
    hcpd_spectrum,
    get_hcpd_spectrum(src = hcpd_spectrum_src, by_run = by_run, hcpd_exclusion)
  ),
  tar_target(
    hcpya_spectrum_src,
    "data/motion/derivatives/hcpya_spectrum.parquet",
    format = "file"
  ),
  tar_target(
    hcpya_spectrum,
    get_hcpya_spectrum(
      src = hcpya_spectrum_src,
      by_run = by_run,
      excluded = hcpya_exclusion
    )
  ),
  tar_target(
    abcd_spectrum_src,
    "data/motion/derivatives/abcd_spectrum.parquet",
    format = "file"
  ),
  tar_target(
    abcd_spectrum,
    get_abcd_spectrum(
      src = abcd_spectrum_src,
      by_run = by_run,
      excluded = abcd_exclusion
    )
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
      by_run = by_run,
      excluded = spacetop_exclusion
    )
  ),
  tar_target(mriqc_src, "data/bold", format = "file"),
  tar_target(
    age_chart,
    write_png(
      make_fig_demographics(by_run, demographics, mriqc_src),
      "figures/age-chart.png",
      width = 6.5,
      height = 4
    ),
    packages = c("patchwork", "ggplot2"),
    format = "file"
  ),
  tar_target(
    fig_by_time_hcpya_file,
    write_png(
      make_fig_by_time_hcpya(by_time, hcpya_events = hcpya_events),
      "figures/by-time-hcpya.png",
      width = 9.5,
      height = 6
    ),
    packages = c("patchwork"),
    format = "file"
  ),
  tar_target(
    hcpd_event_times,
    "data/events_consolidated/hcp_dev_event_times.csv",
    format = "file"
  ),
  tar_target(
    fig_by_time_hcpd,
    write_png(
      make_fig_by_time_hcpd(
        hcpd = hcpd,
        hcpd_event_times = hcpd_event_times,
        by_time = by_time
      ),
      "figures/by-time-hcpd.png",
      width = 6,
      height = 5
    ),
    packages = c("patchwork"),
    format = "file"
  ),
  tar_target(
    hcpa_event_times,
    "data/events_consolidated/hcp_aging_event_times.csv",
    format = "file"
  ),
  tar_target(
    fig_by_time_hcpa,
    write_png(
      make_fig_by_time_hcpa(
        hcpa = hcpa,
        hcpa_event_times = hcpa_event_times,
        by_time = by_time
      ),
      "figures/by-time-hcpa.png",
      width = 6,
      height = 5
    ),
    packages = c("patchwork"),
    format = "file"
  ),
  tar_target(
    fig_by_time_ukb,
    write_png(
      make_fig_by_time_ukb(
        ukb_events = ukb_events,
        ukb_responses = ukb_responses,
        by_time = by_time
      ),
      "figures/by-time-ukb.png",
      width = 5,
      height = 2
    ),
    packages = c("patchwork"),
    format = "file"
  ),
  tar_target(
    abcd_events2,
    "data/events_consolidated/abcd_events.tsv",
    format = "file"
  ),
  tar_target(
    fig_by_time_abcd,
    write_png(
      make_fig_by_time_abcd(
        abcd_events = abcd_events2,
        by_time = by_time
      ),
      "figures/by-time-abcd.png",
      width = 5,
      height = 4
    ),
    packages = c("patchwork"),
    format = "file"
  ),
  tar_target(
    fig_by_time_spacetop,
    write_png(
      make_fig_by_time_spacetop(
        by_time = by_time
      ),
      "figures/by-time-spacetop.png",
      width = 5,
      height = 3
    ),
    packages = c("patchwork"),
    format = "file"
  ),
  tar_target(
    fig_by_run,
    write_png(
      make_fig_by_run(
        by_run = by_run
      ),
      "figures/by-run.png",
      width = 6.5,
      height = 7.7
    ),
    packages = c("patchwork"),
    format = "file"
  ),
  tar_target(
    fig_all_motion_exclusion,
    write_png(
      make_fig_all_motion_exclusion(lost = lost),
      "figures/all_motion_exclusion.png",
      width = 6.5,
      height = 6.5
    ),
    packages = c("patchwork"),
    format = "file"
  ),
  tar_target(
    fig_all_motion_exclusion2,
    write_png(
      make_fig_all_motion_exclusion(lost = lost, filtered = TRUE),
      "figures/all_motion_exclusion_filtered.png",
      width = 6.5,
      height = 6.5
    ),
    packages = c("patchwork"),
    format = "file"
  ),
  tar_target(
    fig_cluster,
    write_png(
      make_fig_cluster(by_run),
      "figures/cluster.png",
      width = 6,
      height = 5
    ),
    packages = c("patchwork"),
    format = "file"
  ),
  tar_target(
    fig_lost_by_group,
    write_png(
      make_fig_lost_by_group(by_run, demographics),
      "figures/lost_by_group.png",
      width = 6,
      height = 5
    ),
    packages = c("patchwork"),
    format = "file"
  )
)
