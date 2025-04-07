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

targets::tar_option_set(
  trust_timestamps = TRUE,
  format = "parquet"
)


list(
  tar_target(ukb_source, "data/motion/derivatives/ukb.parquet", format = "file"),
  tar_target(ukb_exclusion_src, "data/exclusion/ukb_exclusion.tsv", format = "file"),
  tar_target(
    ukb_exclusion,
    get_ukb_exclusion(ukb_exclusion_src),
    deployment = "main"
  ),
  tar_target(ukb, get_ukb(ukb_source, ukb_exclusion)),
  tar_target(ukb_beh, "data/ukb677207_bulk.parquet", format = "file"),
  tar_target(ukb_design_mat, "data/ukb/design.mat", format = "file"),
  # tar_target(ukb_events, get_ukb_design(ukb_design_mat)),
  tar_target(ukb_eprime, "data/1000513_25748_2_0.txt", format = "file"),
  tar_target(ukb_responses, get_ukb_responses(ukb_eprime)),
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
    "data/motion/derivatives/human-connectome-project-openaccess.parquet",
    format = "file"
  ),
  tar_target(hcpya_exclusion, get_hcpya_exclusion(hcpya_source)),
  tar_target(
    hcpya,
    get_hcp(hcpya_source, hcpya_exclusion),
    deployment = "main"
  ),
  # tar_target(hcpya_events, get_hcpya_events("data/hcp_evs")),
  tar_target(
    hcpa_source,
    "data/motion/derivatives/HCPAgingRec.parquet",
    format = "file"
  ),
  tar_target(hcpa, get_hcp(hcpa_source), deployment = "main"),
  tar_target(
    hcpd_source,
    "data/motion/derivatives/HCPDevelopmentRec.parquet",
    format = "file"
  ),
  tar_target(
    hcpd_exclusion_src,
    "data/exclusion/hcp_dev.csv",
    format = "file",
    deployment = "main"
  ),
  tar_target(hcpd_exclusion, get_hcpd_exclusion(hcpd_exclusion_src)),
  tar_target(hcpd, get_hcp(hcpd_source, hcpd_exclusion)),
  tar_target(
    abcd_source,
    "data/motion/derivatives/abcd.parquet",
    format = "file"
  ),
  tar_target(abcd_exclusion, get_abcd_exclusion(abcd_source, abcd_demographics)),
  # tar_target(abcd_events, get_abcd_events()),
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
  tar_target(ukb_demographics, get_ukb_demographics()),
  tar_target(hcpya_demographics, get_hcpya_demographics()),
  tar_target(spacetop_demographics, get_spacetop_demographics()),
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
    hcpa_spectrum,
    get_hcpa_spectrum(
      src = "data/motion/derivatives/hcpa_spectrum.parquet",
      by_run = by_run
    )
  ),
  tar_target(
    hcpd_spectrum,
    get_hcpd_spectrum(
      src = "data/motion/derivatives/hcpd_spectrum.parquet",
      by_run = by_run,
      hcpd_exclusion
    )
  ),
  tar_target(
    hcpya_spectrum,
    get_hcpya_spectrum(
      src = "data/motion/derivatives/hcpya_spectrum.parquet",
      by_run = by_run,
      excluded = hcpya_exclusion
    )
  ),
  tar_target(
    abcd_spectrum,
    get_abcd_spectrum(
      src = "data/motion/derivatives/abcd_spectrum.parquet",
      by_run = by_run,
      excluded = abcd_exclusion
    )
  ),
  tar_target(
    spacetop_spectrum,
    get_spacetop_spectrum(
      src = "data/motion/derivatives/spacetop_spectrum.parquet",
      by_run = by_run,
      excluded = spacetop_exclusion
    )
  ),
  tar_target(
    qc_fd_hcpya,
    get_qc_fd_hcpya(hcpya = hcpya, by_run = by_run, n_iter = 100, n_sub = 50),
    deployment = "main"
  ),
  tarchetypes::tar_group_by(
    qc_fd_hcpya_sub,
    dplyr::filter(qc_fd_hcpya, iter == 0),
    sub,
    deployment = "main"
  ),
  tarchetypes::tar_group_by(
    qc_fd_hcpya_iter,
    dplyr::group_nest(qc_fd_hcpya, iter, sub, filtered),
    iter, sub, filtered,
    deployment = "main"
  ),
  tar_target(cleaned, c(TRUE, FALSE), format = "qs", deployment = "main"),
  tar_target(
    qcs_cor_hcpya_gold,
    get_cor_by_thresh_hcpya_gold(qc_fd_hcpya_sub, cleaned = cleaned),
    pattern = cross(map(qc_fd_hcpya_sub), cleaned)
  )
)
