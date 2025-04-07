get_hcp <- function(src, exclusion = NULL) {
  d <- arrow::open_dataset(src) |>
    dplyr::filter(t > 0) |>
    dplyr::collect() |>
    add_run() |>
    dplyr::mutate(
      time = t * 0.72,
      task = stringr::str_to_lower(task),
      run = dplyr::case_match(
        task,
        "rest1" ~ 1L,
        "rest2" ~ 2L,
        .default = run
      ),
      task = dplyr::case_when(
        task == "rest1a" ~ "resta",
        task == "rest1b" ~ "restb",
        task == "rest2a" ~ "resta",
        task == "rest2b" ~ "restb",
        stringr::str_detect(task, "rest") ~ "rest",
        .default = task
      )
    ) |>
    do_casting()

  if (stringr::str_detect(src, "openaccess")) {
    out <- readr::read_csv(
      fs::dir_ls(here::here("data/sessionSummaryCSV_1200Release"), glob = "*csv"),
      id = "sub",
      show_col_types = FALSE
    ) |>
      dplyr::filter(stringr::str_detect(`Scan Type`, "fMRI$")) |>
      dplyr::filter(stringr::str_detect(`Scan Description`, "MOVIE|RET", TRUE)) |>
      dplyr::mutate(
        sub = stringr::str_extract(sub, "[[:digit:]]{6}"),
        task = stringr::str_extract(`Scan Description`, "WM|REST|LANGUAGE|SOCIAL|GAMBLING|EMOTION|MOTOR|RELATIONAL") |>
          stringr::str_to_lower(),
        run = dplyr::case_when(
          stringr::str_detect(`Scan Description`, "REST1") ~ 1L,
          stringr::str_detect(`Scan Description`, "REST2") ~ 2L,
          TRUE ~ 1L,
        ),
        ped = stringr::str_extract(`Scan Description`, "(?<=_)(RL|LR|AP|PA)"),
        `Session Day` = stringr::str_extract(`Session Day`, "[[:digit:]]")
      ) |>
      dplyr::filter(stringr::str_detect(ped, "A|P", TRUE)) |>
      dplyr::mutate(
        scan = rank(`Acquisition Time`),
        .by = c(sub, task, `Session Day`)
      ) |>
      dplyr::mutate(
        scan = scan + 2 * (as.integer(factor(`Session Day`)) - 1),
        .by = c(sub, task)
      ) |>
      dplyr::select(sub, task, run, ped, ses = `Session Day`, scan) |>
      dplyr::right_join(
        dplyr::select(d, -ses),
        by = dplyr::join_by(sub, task, ped, run)
      ) |>
      dplyr::filter(
        !(sub == 196952 & task == "wm"),
        !(sub == 748662 & task == "social"),
        !(sub == 809252 & task == "social"),
        !(sub %in% c(103010, 113417, 116423, 120010, 121719, 127226, 130114, 143830, 169040, 185038, 189652, 202820, 204218, 329844, 385046, 401422, 462139, 469961, 644246, 688569, 723141, 908860, 908860, 969476, 971160, 168139, 144428, 168139, 186545, 192237, 223929, 320826, 644044, 822244, 870861, 947668))
      ) |>
      dplyr::mutate(ses = "1")
    # https://wiki.humanconnectome.org/docs/HCP%20Data%20Release%20Updates%20Known%20Issues%20and%20Planned%20fixes.html
  } else if (stringr::str_detect(src, "Aging")) {
    out <- get_hcpad_ses("data/hcpa_sessions") |>
      dplyr::right_join(
        dplyr::select(d, -ses),
        by = dplyr::join_by(sub, task, ped, run)
      ) |>
      dplyr::mutate(ses = "1")
  } else {
    # https://www-sciencedirect-com.proxy1.library.jhu.edu/science/article/pii/S1053811918318652?via%3Dihub
    # table 1
    out <- get_hcpad_ses("data/hcpd_sessions") |>
      dplyr::right_join(
        dplyr::select(d, -ses),
        by = dplyr::join_by(sub, task, ped, run)
      ) |>
      dplyr::mutate(ses = "1") |>
      dplyr::filter(!(is.na(scan) & task == "restb")) |> # unclear why these scans do not appear in sessions.csv
      dplyr::filter(!(is.na(scan) & task == "rest" & run == 2))
  }
  if (!is.null(exclusion)) {
    out <- dplyr::anti_join(out, exclusion)
  }
  out
}

get_hcpad_ses <- function(srcs) {
  readr::read_csv(
    fs::dir_ls(here::here(srcs), glob = "*csv", recurse = TRUE),
    id = "sub",
    show_col_types = FALSE
  ) |>
    dplyr::filter(stringr::str_detect(`Scan Type`, "fMRI$")) |>
    dplyr::filter(!(is.na(`Percent Complete`) & stringr::str_detect(`Series Description`, "REST"))) |>
    dplyr::filter(!(is.na(CDB_Description))) |>
    dplyr::mutate(
      sub = stringr::str_extract(sub, "[[:digit:]]{6,7}") |>
        as.integer() |>
        as.character(),
      task = stringr::str_extract(CDB_Description, "REST1a|REST1b|REST2a|REST2b|REST|VISMOTOR|CARIT|FACENAME|EMOTION|GUESSING") |>
        stringr::str_to_lower(),
      task = dplyr::case_when(
        task == "rest1a" ~ "resta",
        task == "rest1b" ~ "restb",
        task == "rest2a" ~ "resta",
        task == "rest2b" ~ "restb",
        stringr::str_detect(task, "rest") ~ "rest",
        .default = task
      ),
      run = dplyr::case_when(
        stringr::str_detect(CDB_Description, "REST1") ~ 1L,
        stringr::str_detect(CDB_Description, "REST2") ~ 2L,
        .default = 1L
      ),
      ped = stringr::str_extract(`Series Description`, "(?<=_)(RL|LR|AP|PA)"),
      `Session Day` = stringr::str_extract(`Session Day`, "[[:digit:]]")
    ) |>
    dplyr::mutate(
      scan = rank(`Acquisition Time`),
      .by = c(sub, task, `Session Day`)
    ) |>
    dplyr::mutate(
      scan = scan + 2 * (as.integer(factor(`Session Day`)) - 1),
      .by = c(sub, task)
    ) |>
    dplyr::mutate(
      scan = dplyr::if_else(task == "restb" & scan == 3, 2, scan)
    ) |>
    dplyr::select(sub, task, run, ped, ses = `Session Day`, scan)
}



get_hcp_design0 <- function(src) {
  .evs <- fs::dir_ls(src, recurse = TRUE, glob = "*.txt")
  .evs <- .evs[stringr::str_detect(.evs, "Sync", TRUE)]
  lengths <- purrr::map_dbl(.evs, fpeek::peek_count_lines)
  .evs <- .evs[lengths > 0]

  batches <- tibble::tibble(
    src = .evs
  ) |>
    dplyr::mutate(
      batch = 1:dplyr::n() %% 10
    ) |>
    dplyr::group_nest(batch)

  purrr::map(
    batches$data,
    ~ readr::read_tsv(
      .x$src,
      col_names = c("event", "onset", "duration", "amplitude"),
      col_select = c("event", "onset", "duration"),
      col_types = "iddd",
      id = "path"
    )
  ) |>
    dplyr::bind_rows() |>
    dplyr::mutate(
      task = stringr::str_extract(
        path,
        "MOTOR|SOCIAL|EMOTION|WM|RELATIONAL|GAMBLING|LANGUAGE"
      ) |>
        factor(
          levels = c("WM", "GAMBLING", "MOTOR", "LANGUAGE", "SOCIAL", "RELATIONAL", "EMOTION"),
          ordered = TRUE
        ),
      ped = stringr::str_extract(path, "RL|LR") |> as.factor(),
      type = stringr::str_extract(path, "([[:alpha:]]+)(?=.txt)")
    ) |>
    dplyr::mutate(sub = stringr::str_extract(path, "[[:digit:]]{6}"))
}

get_hcpya_events <- function(src) {
  hcp_design <- get_hcp_design0(src)

  wm_hcp_design <- hcp_design |>
    dplyr::filter(task == "WM") |>
    dplyr::filter(type %in% c("tools", "places", "faces", "body")) |>
    dplyr::arrange(onset) |>
    dplyr::mutate(
      event = 1:dplyr::n(),
      .by = c(ped, task, sub)
    ) |>
    dplyr::summarise(
      onset = mean(onset),
      duration = mean(duration),
      .by = c(event, task, ped)
    ) |>
    dplyr::mutate(type = "block")

  gambling_hcp_design <- hcp_design |>
    dplyr::filter(task == "GAMBLING") |>
    dplyr::filter(type %in% c("win", "loss")) |>
    dplyr::arrange(onset) |>
    dplyr::mutate(
      event = 1:dplyr::n(),
      .by = c(ped, task, sub)
    ) |>
    dplyr::summarise(
      onset = mean(onset),
      duration = mean(duration),
      .by = c(event, task, ped)
    ) |>
    dplyr::mutate(type = "block")

  social_hcp_stim <- hcp_design |>
    dplyr::filter(task == "SOCIAL") |>
    dplyr::filter(type %in% c("rnd", "mental")) |>
    dplyr::arrange(onset) |>
    dplyr::mutate(
      event = 1:dplyr::n(),
      .by = c(ped, task, sub)
    ) |>
    dplyr::summarise(
      onset = mean(onset),
      duration = mean(duration),
      .by = c(event, task, ped)
    ) |>
    dplyr::mutate(type = "block")

  rel_hcp_stim <- hcp_design |>
    dplyr::filter(task == "RELATIONAL") |>
    dplyr::filter(type %in% c("relation", "match")) |>
    dplyr::arrange(onset) |>
    dplyr::mutate(
      event = 1:dplyr::n(),
      .by = c(ped, task, sub)
    ) |>
    dplyr::summarise(
      onset = mean(onset),
      duration = mean(duration),
      .by = c(event, task, ped)
    ) |>
    dplyr::mutate(type = "block")

  lang_hcp_stim <- hcp_design |>
    dplyr::filter(task == "LANGUAGE") |>
    dplyr::filter(
      type %in% c("story", "math"),
      stringr::str_detect(path, "response", TRUE),
      stringr::str_detect(path, "question", TRUE),
      stringr::str_detect(path, "present", TRUE)
    ) |>
    dplyr::arrange(onset) |>
    dplyr::mutate(
      event = 1:dplyr::n(),
      .by = c(ped, task, sub, type)
    ) |>
    dplyr::summarise(
      onset = mean(onset),
      duration = mean(duration),
      .by = c(event, task, ped, type)
    )


  hcp_design |>
    dplyr::filter(!task %in% c("WM", "GAMBLING", "SOCIAL", "LANGUAGE", "RELATIONAL")) |>
    dplyr::mutate(
      type = dplyr::case_when(
        type == "neut" ~ "shapes",
        type == "fear" ~ "faces",
        TRUE ~ type
      )
    ) |>
    na.omit() |>
    dplyr::summarise(
      onset = mean(onset),
      duration = mean(duration),
      .by = c(event, task, ped, type)
    ) |>
    dplyr::bind_rows(wm_hcp_design, gambling_hcp_design, social_hcp_stim, rel_hcp_stim, lang_hcp_stim)
}

get_exclude_hcp <- function(path) {
  readr::read_csv(
    path,
    col_select = c("sub" = "Subject", "QC_Issue")
  ) |>
    dplyr::filter(!is.na(QC_Issue)) |>
    dplyr::select(sub)
}

get_hcpya_demographics <- function() {
  hcpya_un <- readr::read_csv(here::here("data/unrestricted_martin_2_5_2024_10_18_12.csv")) |>
    dplyr::select(sub = Subject, gender = Gender)

  readr::read_csv(here::here("data/RESTRICTED_martin_2_5_2024_10_18_28.csv")) |>
    dplyr::select(sub = Subject, age = Age_in_Yrs, race = Race, ethnicity = Ethnicity, bmi = BMI) |>
    dplyr::left_join(hcpya_un, by = dplyr::join_by(sub)) |>
    dplyr::mutate(sub = as.character(sub))
}

get_bmi_from_vitals <- function(src) {
  read_nda(src) |>
    dplyr::mutate(
      weight_kg = weight_std * 0.4535924,
      height_m = vtl007 * 0.0254,
      bmi = weight_kg / height_m^2
    ) |>
    dplyr::select(sub = src_subject_id, bmi)
}


get_hcp_aging_demographics <- function() {
  vitals <- get_bmi_from_vitals("data/demographics/HCPAgingRec_vitals01.txt")

  read_nda(here::here("data/demographics/HCPAgingRec_ndar_subject01.txt")) |>
    dplyr::select(
      sub = src_subject_id,
      age = interview_age,
      race = race,
      ethnicity = ethnic_group,
      site = site,
      sex = sex
    ) |>
    dplyr::left_join(vitals) |>
    dplyr::mutate(
      age = age / 12,
      sub = stringr::str_remove(sub, "HCA") |> as.integer() |> as.character()
    )
}

get_hcp_dev_demographics <- function() {
  vitals <- get_bmi_from_vitals("data/demographics/HCPDevelopmentRec_vitals01.txt")
  read_nda(here::here("data/demographics/HCPDevelopmentRec_ndar_subject01.txt")) |>
    dplyr::select(
      sub = src_subject_id,
      age = interview_age,
      race = race,
      ethnicity = ethnic_group,
      site = site,
      sex = sex
    ) |>
    dplyr::left_join(vitals) |>
    dplyr::mutate(
      age = age / 12,
      sub = stringr::str_remove(sub, "HCD") |> as.integer() |> as.character()
    )
}

get_hcpya_exclusion <- function(src) {
  expected <- arrow::open_dataset(src) |>
    dplyr::summarise(n_tr = max(t), .by = c(sub, task, ped)) |>
    dplyr::count(n_tr, task, ped) |>
    dplyr::collect() |>
    dplyr::slice_max(order_by = n, n = 1, with_ties = FALSE, by = c(task, ped)) |>
    dplyr::select(-n)

  arrow::open_dataset(src) |>
    dplyr::summarise(n_tr = max(t), .by = c(sub, task, ped)) |>
    dplyr::collect() |>
    dplyr::anti_join(expected) |>
    dplyr::select(sub, task, ped)
}

get_hcpd_exclusion <- function(src = "data/exclusion/hcp_dev.csv") {
  readr::read_csv(src, col_types = "cccic") |>
    dplyr::mutate(
      task = stringr::str_to_lower(task),
      sub = stringr::str_remove(sub, "HCD") |>
        as.integer() |> as.character()
    )
}
