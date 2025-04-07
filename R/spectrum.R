
rescale <- function(.data, ...){
  .data |>
    dplyr::mutate(
      avg = factor(avg, ordered = TRUE),
      pxx = 10*log10(pxx)) |>
    na.omit() |>
    dplyr::group_nest(..., avg) |>
    dplyr::mutate(
      data = purrr::map(
        data,
        ~tibble::tibble(
          pxx = lm(scale(pxx) ~ pxx, data=.x) |> predict(),
          freq=.x$freq,
        )
      ),
      lower = purrr::map_dbl(data, ~quantile(.x$pxx, 0.01)),
      upper = purrr::map_dbl(data, ~quantile(.x$pxx, 0.96))
    ) |>
    tidyr::unnest(data) |>
    dplyr::mutate(
      pxx = (pxx - min(pxx)) / (max(pxx) - min(pxx)),  
      pxx=dplyr::if_else(pxx<lower, lower, pxx),
      pxx=dplyr::if_else(pxx>upper, upper, pxx), 
      .by = c(...)
    ) |>
    dplyr::select(-lower, -upper)
}

.avg_by_run <- function(by_run){
  by_run |>
    dplyr::filter(!filtered) |>
    dplyr::summarise(
      avg = median(loc),
      .by = c(sub, dataset)
    )
}

.clean <- function(.data, N=100){
  .data |>
    semi_join(count(.data, freq) |> filter(n>N))
}


get_hcpa_spectrum <- function(src, by_run){
  arrow::read_parquet(src) |>
    dplyr::left_join(
      .avg_by_run(dplyr::filter(by_run, dataset=="hcpa")), 
      by = dplyr::join_by(sub)) |>
    dplyr::filter(!(task=="REST1" & ped=="PA")) |>
    rescale(task, ped, sub, param)
}

get_hcpd_spectrum <- function(src, by_run, excluded){
  arrow::read_parquet(src) |>
    dplyr::left_join(
      .avg_by_run(dplyr::filter(by_run, dataset=="hcpd")), 
      by = dplyr::join_by(sub)) |>
    dplyr::filter(!(task%in%c("REST1a","REST1b","REST2a","REST2b"))) |>
    dplyr::mutate(
      task = stringr::str_to_lower(task)
    ) |>
    dplyr::filter(!is.na(avg)) |>
    dplyr::anti_join(excluded) |>    
    rescale(task, ped, sub, param)
}

get_hcpya_spectrum <- function(src, by_run, excluded){
  arrow::read_parquet(src) |>
    dplyr::left_join(
      .avg_by_run(dplyr::filter(by_run, dataset=="hcpya")), 
      by = dplyr::join_by(sub)
    ) |>
    dplyr::filter(!is.na(avg)) |>
    dplyr::anti_join(excluded) |>
    rescale(task, ped, sub, param)
}

get_abcd_spectrum <- function(src, by_run, excluded){
  
  by_run_abcd <- by_run  |>
    dplyr::filter(dataset=="abcd", task=="rest",  !filtered) |>
    dplyr::rename(avg=loc) |> 
    dplyr::anti_join(
      readr::read_table("data/abcd_exclusion.tsv", show_col_types = FALSE),
      by = dplyr::join_by(task, ses, sub)) |>
    dplyr::summarise(avg=mean(avg), .by = sub) 
  arrow::read_parquet(src) |>
    dplyr::anti_join(excluded, by = dplyr::join_by(sub, ses, task, run)) |>
    dplyr::left_join(by_run_abcd) |>
    dplyr::filter(!is.na(avg)) |>
    rescale(task, sub, ses, run, param) |>
    dplyr::mutate(
      ses = factor(
        ses,
        levels = c(
          "baselineYear1Arm1", 
          "2YearFollowUpYArm1", 
          "4YearFollowUpYArm1"
        ),
        ordered = TRUE
      ),
      ses=forcats::fct_recode(
        ses,
        baseline="baselineYear1Arm1",
        Year2="2YearFollowUpYArm1",
        Year4="4YearFollowUpYArm1"
      ) |>
        forcats::fct_relevel("Year2", after=Inf) |>
        forcats::fct_relevel("Year4", after=Inf)
    )
  
}

get_spacetop_spectrum <- function(src, by_run, excluded){
  arrow::read_parquet(src) |>
    dplyr::left_join(
      .avg_by_run(dplyr::filter(by_run, dataset=="spacetop")), 
      by = dplyr::join_by(sub)
    ) |>
    dplyr::filter(!is.na(avg)) |>
    dplyr::anti_join(excluded, by=dplyr::join_by(sub, ses, task, run)) |>
    rescale(task, ses, run, sub, param)
}

