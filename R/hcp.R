
get_hcp <- function(sources, exclude_hcp){
  batches <- tibble::tibble(
    src = sources
  ) |>
    dplyr::mutate(
      batch = 1:dplyr::n() %% 10
    ) |>
    dplyr::group_nest(batch)
  
  purrr::map(
    batches$data,
    ~readr::read_tsv(
      .x$src,
      col_names = c(
        "t",
        "trans_x", "trans_y", "trans_z",
        "rot_x", "rot_y", "rot_z",
        "trans_x_derivative1", 
        "trans_y_derivative1", 
        "trans_z_derivative1",
        "rot_x_derivative1", 
        "rot_y_derivative1", 
        "rot_z_derivative1"),
      col_types = "idddddddddddd",
      id = "path")) |>
    dplyr::bind_rows() |>
    dplyr::mutate(
      dplyr::across(tidyselect::starts_with("rot"), deg_2_rad),
      sub = stringr::str_extract(path, "[[:digit:]]{6}") |>
        as.integer(),
      ped = stringr::str_extract(path, "LR|RL") |>
        as.factor(),
      task = stringr::str_extract(
        path, 
        "MOTOR|SOCIAL|EMOTION|WM|RELATIONAL|GAMBLING|LANGUAGE") |>
        factor(
          levels = c("WM", "GAMBLING", "MOTOR", "LANGUAGE", "SOCIAL", "RELATIONAL", "EMOTION"),
          ordered = TRUE
        ),
      ses = dplyr::case_match(
        task,
        c("WM", "GAMBLING", "MOTOR") ~ 1, # SB before two REST scans
        c("LANGUAGE", "SOCIAL", "RELATIONAL", "EMOTION") ~ 2, # SB Immediately before
      ),
      .by = c(path)
    ) |>
    dplyr::select(-path) |>
    dplyr::anti_join(exclude_hcp, by = dplyr::join_by(sub))
}


get_hcp_design0 <- function(src){
  
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
    ~readr::read_tsv(
      .x$src,
      col_names = c("event", "onset", "duration", "amplitude"),
      col_select = c("event","onset", "duration"),
      col_types = "iddd",
      id = "path")) |>
    dplyr::bind_rows() |>
    dplyr::mutate(
      task = stringr::str_extract(
        path,
        "MOTOR|SOCIAL|EMOTION|WM|RELATIONAL|GAMBLING|LANGUAGE")|>
        factor(
          levels = c("WM", "GAMBLING", "MOTOR", "LANGUAGE", "SOCIAL", "RELATIONAL", "EMOTION"),
          ordered = TRUE
        ),
      ped = stringr::str_extract(path, "RL|LR") |> as.factor(),
      type = stringr::str_extract(path, "([[:alpha:]]+)(?=.txt)")
    ) |>
    dplyr::mutate(sub = stringr::str_extract(path, "[[:digit:]]{6}")) 
    
}

get_hcp_design <- function(hcp_design){
  
  wm_hcp_design <- hcp_design  |>
    dplyr::filter(task=="WM") |>
    dplyr::filter(type %in% c("tools", "places", "faces", "body")) |>
    dplyr::arrange(onset) |>
    dplyr::mutate(
      event = 1:dplyr::n(),
      .by = c(ped, task, sub)
    ) |>
    dplyr::summarise(
      onset = mean(onset / 0.72),
      duration = mean(duration / 0.72),
      .by = c(event, task, ped)
    ) |>
    dplyr::mutate(type = "block")
  
  gambling_hcp_design <- hcp_design  |>
    dplyr::filter(task=="GAMBLING") |>
    dplyr::filter(type %in% c("win", "loss")) |>
    dplyr::arrange(onset) |>
    dplyr::mutate(
      event = 1:dplyr::n(),
      .by = c(ped, task, sub)
    ) |>
    dplyr::summarise(
      onset = mean(onset / 0.72),
      duration = mean(duration / 0.72),
      .by = c(event, task, ped)
    ) |>
    dplyr::mutate(type = "block")
  
  social_hcp_stim <- hcp_design  |>
    dplyr::filter(task=="SOCIAL") |>
    dplyr::filter(type %in% c("rnd", "mental")) |>
    dplyr::arrange(onset) |>
    dplyr::mutate(
      event = 1:dplyr::n(),
      .by = c(ped, task, sub)
    ) |>
    dplyr::summarise(
      onset = mean(onset / 0.72),
      duration = mean(duration / 0.72),
      .by = c(event, task, ped)
    ) |>
    dplyr::mutate(type = "block")
  
  rel_hcp_stim <- hcp_design  |>
    dplyr::filter(task=="RELATIONAL") |>
    dplyr::filter(type %in% c("relation", "match")) |>
    dplyr::arrange(onset) |>
    dplyr::mutate(
      event = 1:dplyr::n(),
      .by = c(ped, task, sub)
    ) |>
    dplyr::summarise(
      onset = mean(onset / 0.72),
      duration = mean(duration / 0.72),
      .by = c(event, task, ped)
    ) |>
    dplyr::mutate(type = "block")
  
  lang_hcp_stim <- hcp_design  |>
    dplyr::filter(task=="LANGUAGE") |>
    dplyr::filter(
      type %in% c("story", "math"), 
      stringr::str_detect(path, "response", TRUE),
      stringr::str_detect(path, "question", TRUE),
      stringr::str_detect(path, "present", TRUE)) |>
    dplyr::arrange(onset) |>
    dplyr::mutate(
      event = 1:dplyr::n(),
      .by = c(ped, task, sub, type)
    ) |>
    dplyr::summarise(
      onset = mean(onset / 0.72),
      duration = mean(duration / 0.72),
      .by = c(event, task, ped, type)
    ) 
  
  
  hcp_design |>
    dplyr::filter(!task%in%c("WM" ,"GAMBLING" ,"SOCIAL", "LANGUAGE",  "RELATIONAL")) |>
    dplyr::mutate(
      type = dplyr::case_when(
        type == "neut" ~ "shapes",
        type == "fear" ~ "faces",
        TRUE ~ type
      )
    ) |>
    na.omit() |>
    dplyr::summarise(
      onset = mean(onset / 0.72),
      duration = mean(duration / 0.72),
      .by = c(event, task, ped, type)
    ) |>
    dplyr::bind_rows(wm_hcp_design, gambling_hcp_design, social_hcp_stim,rel_hcp_stim, lang_hcp_stim) 
  
}

get_exclude_hcp <- function(path){
  readr::read_csv(
    path,
    col_select = c("sub"="Subject", "QC_Issue")) |>
    dplyr::filter(!is.na(QC_Issue)) |>
    dplyr::select(sub)
}

