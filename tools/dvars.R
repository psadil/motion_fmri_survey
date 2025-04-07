library(tidyr)
library(dplyr)
library(ggplot2)
library(purrr)
library(stringr)

hcpya_not_avail <- c(
  110613, 113417, 113821, 120010, 121719, 130518, 139637, 143830, 146836,
  168139, 175035, 176239, 185038, 189652, 199958, 201515, 202820, 385046,
  401422, 415837, 433839, 462139, 465852, 469961, 644246, 656657, 688569,
  723141, 767464, 872764, 943862, 965367, 969476, 987983, 994273, 433839,
  # https://wiki.humanconnectome.org/docs/HCP%20Data%20Release%20Updates%20Known%20Issues%20and%20Planned%20fixes.html
  # 3T Functional Preprocessing Error of all 3T “RL” fMRI runs in 25 Subjects
  103010, 113417, 116423, 120010, 121719, 127226, 130114, 143830, 169040,
  185038, 189652, 202820, 204218, 329844, 385046, 401422, 462139, 469961,
  644246, 688569, 723141, 908860, 943862, 969476, 971160,
  196952, 748662, 809252, 144428, 186545, 192237, 223929, 320826, 644044,
  822244, 870861, 947668,
  # Subjects without Field Maps for Structural scans
  102614, 111009, 111514, 115017, 121416, 130821, 138332, 179952, 299760,
  300618, 392750, 406432, 429040, 633847, 662551, 679770, 688569, 693461,
  815247,
  # duplicated name
  142626
)



ukb_dvars <- arrow::open_dataset("data/dvars/dataset=ukb", format = "arrow") |>
  mutate(
    keep = (stringr::str_detect(src, "20249") & t < 332) |
      (stringr::str_detect(src, "20227") & t < 490)
  ) |>
  filter(ses == 2, t > 0, keep, !src == "20227-filtered_func_data") |>
  summarise(
    DPD = mean(DPD),
    ZD = mean(ZD),
    .by = c(t, src)
  ) |>
  collect() |>
  mutate(
    datatype = str_extract(src, "[[:digit:]]{5}")
  ) |>
  select(-src)



p <- ukb_dvars |>
  pivot_longer(-c(t, datatype)) |>
  ggplot(aes(x = t, y = value)) +
  facet_grid(name ~ datatype, scales = "free") +
  geom_line()

ggsave("figures/ukb_dvars.png", p, width = 6)

hcpya_dvars <- arrow::open_dataset("data/dvars/dataset=human-connectome-project-openaccess/", format = "arrow") |>
  filter(str_detect(src, "MSM")) |>
  filter(!sub %in% hcpya_not_avail) |>
  summarise(
    DPD = mean(DPD),
    ZD = mean(ZD),
    .by = c(t, src)
  ) |>
  collect() |>
  mutate(
    ped = str_extract(src, "RL|LR"),
    task = str_extract(src, "WM|REST1|REST2|EMOTION|RELATIONAL|MOTOR|GAMBLING|SOCIAL|LANGUAGE"),
    ica_cleaned = str_detect(src, "clean")
  )



p <- hcpya_dvars |>
  ggplot(aes(x = t, y = DPD, color = ped, linetype = ica_cleaned)) +
  facet_wrap(~task, scales = "free_x") +
  geom_line(alpha = 0.6)

ggsave("figures/hcpya_dpd.png", p, width = 6)

p <- hcpya_dvars |>
  ggplot(aes(x = t, y = ZD, color = ped, linetype = ica_cleaned)) +
  facet_wrap(~task, scales = "free_x") +
  geom_line(alpha = 0.6)

ggsave("figures/hcpya_zd.png", p, width = 6)


hcpya_sub <- arrow::open_dataset("data/dvars/dataset=human-connectome-project-openaccess/", format = "arrow") |>
  filter(sub == 250427) |>
  filter(str_detect(src, "MSM"), str_detect(src, "REST1")) |>
  select(t, src, sub, DPD, D, value, A) |>
  collect() |>
  mutate(
    ped = str_extract(src, "RL|LR"),
    task = str_extract(src, "WM|REST1|REST2|EMOTION|RELATIONAL|MOTOR|GAMBLING|SOCIAL|LANGUAGE"),
    ica_cleaned = str_detect(src, "clean")
  )




arrow::open_dataset(
  "data/dvars/dataset=human-connectome-project-openaccess/",
  format = "arrow"
) |>
  filter(sub == 250427) |>
  filter(str_detect(src, "MSM")) |>
  select(t, src, sub, DPD) |>
  collect() |>
  mutate(
    ped = str_extract(src, "RL|LR"),
    task = str_extract(src, "WM|REST1|REST2|EMOTION|RELATIONAL|MOTOR|GAMBLING|SOCIAL|LANGUAGE"),
    ica_cleaned = str_detect(src, "clean")
  )



hcpaging_dvars <- arrow::open_dataset("data/dvars/dataset=HCPAgingRec", format = "arrow") |>
  filter(str_detect(src, "Atlas_hp_preclean", negate = TRUE)) |>
  summarise(
    DPD = mean(DPD),
    ZD = mean(ZD),
    .by = c(t, src)
  ) |>
  collect() |>
  mutate(
    ped = str_extract(src, "AP|PA"),
    task = str_extract(src, "REST1|REST2|VISMOTOR|FACENAME|CARIT")
  )

p <- hcpaging_dvars |>
  ggplot(aes(x = t, y = DPD, color = ped)) +
  facet_wrap(~task, scales = "free_x") +
  geom_line(alpha = 0.6)

ggsave("figures/HCPAgingRec_dpd.png", p, width = 6)

p <- hcpaging_dvars |>
  ggplot(aes(x = t, y = ZD, color = ped)) +
  facet_wrap(~task, scales = "free_x") +
  geom_line(alpha = 0.6)

ggsave("figures/HCPAgingRec_zd.png", p, width = 6)


hcpdev_dvars <- arrow::open_dataset("data/dvars/dataset=HCPDevelopmentRec", format = "arrow") |>
  filter(str_detect(src, "Atlas_hp_preclean", negate = TRUE)) |>
  filter(str_detect(src, "REST1a|REST1b|REST2a|REST2b", negate = TRUE)) |>
  summarise(
    DPD = mean(DPD),
    ZD = mean(ZD),
    .by = c(t, src)
  ) |>
  collect() |>
  mutate(
    ped = str_extract(src, "AP|PA"),
    task = str_extract(src, "REST1a|REST1b|REST2a|REST2b|REST1|REST2|EMOTION|GUESSING|CARIT")
  )

p <- hcpdev_dvars |>
  ggplot(aes(x = t, y = DPD, color = ped)) +
  facet_wrap(~task, scales = "free_x") +
  geom_line(alpha = 0.6)

ggsave("figures/HCPDevelopmentRec_dpd.png", p, width = 6)

p <- hcpdev_dvars |>
  ggplot(aes(x = t, y = ZD, color = ped)) +
  facet_wrap(~task, scales = "free_x") +
  geom_line(alpha = 0.6)

ggsave("figures/HCPDevelopmentRec_zd.png", p, width = 6)




abcd_dvars <- arrow::open_dataset("/Users/psadil/data/dvars/derivatives/dvars/dataset=abcd", format = "arrow") |>
  mutate(src = stringr::str_sub(src, start = 21)) |>
  summarise(
    DPD = mean(DPD),
    ZD = mean(ZD),
    .by = c(t, ses, src)
  ) |>
  collect() |>
  mutate(
    task = str_extract(src, "(?<=task-)[[:alnum:]]+(?=_run)"),
    run = str_extract(src, "(?<=run-)[[:alnum:]]+(?=_bold)"),
    .by = c(src),
  ) |>
  mutate(
    ses = factor(
      ses,
      levels = c("baselineYear1Arm1", "2YearFollowUpYArm1", "4YearFollowUpYArm1"),
      ordered = TRUE
    )
  ) |>
  filter(
    t > 10,
    !run %in% c("05", "06"),
    str_detect(task, "rest", TRUE) | (str_detect(task, "rest", FALSE) & (t <= 380))
  )

p <- abcd_dvars |>
  ggplot(aes(x = t, y = DPD, color = run)) +
  facet_grid(task ~ ses, scales = "free_x") +
  geom_line(alpha = 0.6)

ggsave("figures/abcd_dpd.png", p, width = 6)

p <- abcd_dvars |>
  ggplot(aes(x = t, y = ZD, color = run)) +
  facet_grid(task ~ ses, scales = "free_x") +
  geom_line(alpha = 0.6)

ggsave("figures/abcd_zd.png", p, width = 6)


# spacetop


spacetop_dvars <- arrow::open_dataset("data/dvars/dataset=spacetop", format = "arrow") |>
  mutate(src = stringr::str_sub(src, start = 21)) |>
  summarise(
    DPD = mean(DPD),
    ZD = mean(ZD),
    .by = c(t, ses, src)
  ) |>
  collect() |>
  mutate(
    task = str_extract(src, "(?<=task-)[[:alnum:]]+(?=_run)"),
    run = str_extract(src, "(?<=run-)[[:alnum:]]+(?=_bold)"),
    .by = c(src),
  )

p <- spacetop_dvars |>
  ggplot(aes(x = t, y = DPD, color = run)) +
  facet_grid(task ~ ses, scales = "free_x") +
  geom_line(alpha = 0.6)

ggsave("figures/spacetop_dpd.png", p, width = 6)

p <- spacetop_dvars |>
  ggplot(aes(x = t, y = ZD, color = run)) +
  facet_grid(task ~ ses, scales = "free_x") +
  geom_line(alpha = 0.6)

ggsave("figures/spacetop_zd.png", p, width = 6)
