library(dplyr)

d <- nio::to_tbl("~/Desktop/sub-2820279_ses-3_task-rest_bold.nii.gz")

rot_mat <- function(rot_x, rot_y, rot_z) {
  R_x <- matrix(c(
    1, 0, 0, 0,
    0, cos(rot_x), sin(rot_x), 0,
    0, -sin(rot_x), cos(rot_x), 0,
    0, 0, 0, 1
  ), nrow = 4)

  R_y <- matrix(c(
    cos(rot_y), 0, -sin(rot_y), 0,
    0, 1, 0, 0,
    sin(rot_y), 0, cos(rot_y), 0,
    0, 0, 0, 1
  ), nrow = 4)

  R_z <- matrix(c(
    cos(rot_z), sin(rot_z), 0, 0,
    -sin(rot_z), cos(rot_z), 0, 0,
    0, 0, 1, 0,
    0, 0, 0, 1
  ), nrow = 4)

  R_x %*% R_y %*% R_z
}

rot_mat(0, -0, 0.000124956)


n <- RNifti::niftiHeader("~/Desktop/rmsd/sub-2820279_ses-3_task-rest_bold.nii.gz")
aff <- rbind(n$srow_x, n$srow_y, n$srow_z)
centre <- matrix(c(0, 0, 0), ncol = 1)
aff[1:3, 1:3] %*% centre + aff[1:3, 4] - centre


d <- arrow::open_dataset("~/data/motion/derivatives/ukb.parquet") |>
  filter(ses == "3") |>
  filter(sub == "2820279") |>
  filter(task == "rest") |>
  select(-ends_with("filtered")) |>
  collect()

d <- readr::read_delim(
  "/Users/psadil/Desktop/rmsd/sub-2820279_ses-3_task-rest_bold_mcf.par",
  delim = "  ",
  col_names = c("rot_x", "rot_y", "rot_z", "trans_x", "trans_y", "trans_z")
) |>
  select(-X7)
