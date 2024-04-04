
library(targets)
source("R/hcp.R")
source("R/ukb.R")
source("R/mriqc.R")
source("R/utils.R")


list(
  tar_target(
    ukb_sources, 
    fs::dir_ls("data/ukb", recurse = TRUE, glob="*par"), 
    format = "file_fast"),
  tar_target(ukb, get_ukb(ukb_sources), format = "parquet"),
  tar_target(ukb_beh, "data/ukb677207_bulk.parquet", format = "file"),
  tar_target(include_ukb, get_include_ukb(ukb_beh), format = "parquet"),
  tar_target(ukb_reg, get_ukb_reg(ukb, include_ukb), format = "parquet"),
  tar_target(ukb_design_mat, "data/ukb/design.mat", format = "file"),
  tar_target(ukb_design, get_ukb_design(ukb_design_mat), format = "parquet"),
  tar_target(ukb_eprime, "data/1000513_25748_2_0.txt", format = "file"),
  tar_target(ukb_responses, get_ukb_responses(ukb_eprime), format = "parquet"),
  tar_target(
    hcp_sources, 
    fs::dir_ls("data/hcp", recurse = TRUE, glob="*txt"), 
    format = "file_fast"),
  tar_target(hcp_beh, "data/hcp_beh_2_8_2024_15_20_20.csv", format = "file"),
  tar_target(exclude_hcp, get_exclude_hcp(hcp_beh), format = "parquet"),
  tar_target(hcp, get_hcp(hcp_sources, exclude_hcp), format = "parquet"),
  tar_target(hcp_design0, get_hcp_design0("data/hcp_evs"), format = "parquet"),
  tar_target(hcp_design, get_hcp_design(hcp_design0), format = "parquet"),
  tar_target(mriqc_bold, get_mriqc_bold(), format = "parquet")
)
