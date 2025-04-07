library(dplyr)
library(rvest)
library(tidyr)

fields_surf <- xml2::read_html("https://biobank.ndph.ox.ac.uk/showcase/label.cgi?id=197") |>
  html_elements(".listing") |>
  html_text2() |>
  readr::read_tsv() |>
  na.omit() |>
  mutate(
    fields2 = glue::glue("f.{`Field ID`}.2.0"),
    fields3 = glue::glue("f.{`Field ID`}.3.0")
  )

fields_dti <- xml2::read_html("https://biobank.ndph.ox.ac.uk/showcase/label.cgi?id=135") |>
  html_elements(".listing") |>
  html_text2() |>
  readr::read_tsv() |>
  na.omit() |>
  mutate(
    fields2 = glue::glue("f.{`Field ID`}.2.0"),
    fields3 = glue::glue("f.{`Field ID`}.3.0")
  )


d <- arrow::open_dataset("data/ukb677207_bulk.parquet") |>
  select(
    sub = eid,
    any_of(fields_surf$fields2),
    any_of(fields_surf$fields3),
    any_of(fields_dti$fields2),
    any_of(fields_dti$fields3)
  ) |>
  collect() |>
  pivot_longer(
    starts_with("f"),
    names_to = c("f", "field", "visit", "array"),
    names_sep = "\\."
  )


out <- d |>
  select(-f, -array) |>
  left_join(fields_surf, by = join_by(field == `Field ID`)) |>
  left_join(fields_dti, by = join_by(field == `Field ID`), suffix = c(".surf", ".dti")) |>
  mutate(Description = if_else(is.na(Description.surf), Description.dti, Description.surf))


out |>
  select(sub, field, visit, value, Description) |>
  readr::write_tsv("ukb.tsv")
