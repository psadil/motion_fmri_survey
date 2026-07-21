#!/usr/bin/env Rscript
# Generate the bioRxiv author submission sheet from _authors.yml.
#
#   Rscript R/make_author_tsv.R
#
# _authors.yml is the single source of truth: it is also what the manuscript
# itself reads (via `metadata-files` in _quarto.yml), so the submission sheet
# and the rendered PDF can never drift apart.
#
# Column names and their order are taken from author_template.tsv, so if
# bioRxiv changes the template, drop in the new one and re-run.

SRC <- "_authors.yml"
TEMPLATE <- "tools/author_template.tsv"
OUT <- "authors_biorxiv.tsv"

`%||%` <- function(x, y) if (is.null(x)) y else x

meta <- yaml::read_yaml(SRC)

# Affiliation id -> display name. Names are written for LaTeX, so "\&" has to
# be turned back into a plain ampersand for the sheet.
aff_lookup <- vapply(
  meta$affiliations,
  function(a) gsub("\\\\&", "&", a$name),
  character(1)
)
names(aff_lookup) <- vapply(meta$affiliations, function(a) a$id, character(1))

# Split "Ann S. Choe" into first / middle / last.
#
# Family name is the final whitespace-separated token; everything before it is
# the given name. Within the given name, trailing initials ("S.", "V.") are
# treated as middle names, so compound given names such as "Mary Beth" stay
# intact in the first-name column. Any author can override the guess by setting
# first-name / middle-name / last-name in _authors.yml.
split_name <- function(author) {
  parts <- strsplit(trimws(author$name), "\\s+")[[1]]
  last <- parts[length(parts)]
  given <- parts[-length(parts)]

  is_initial <- grepl("^[[:upper:]]\\.$", given)
  # Only trailing initials count as middle names.
  n_mid <- 0L
  while (n_mid < length(given) && is_initial[length(given) - n_mid]) {
    n_mid <- n_mid + 1L
  }
  middle <- if (n_mid > 0L) {
    given[seq(length(given) - n_mid + 1L, length(given))]
  } else {
    character(0)
  }
  first <- if (n_mid > 0L) given[seq_len(length(given) - n_mid)] else given

  list(
    first = author[["first-name"]] %||% paste(first, collapse = " "),
    middle = author[["middle-name"]] %||% paste(middle, collapse = " "),
    last = author[["last-name"]] %||% last
  )
}

rows <- lapply(meta$author, function(a) {
  nm <- split_name(a)
  refs <- vapply(a$affiliations, function(x) x$ref, character(1))
  data.frame(
    Email = a$email %||% "",
    Institution = paste(aff_lookup[refs], collapse = "; "),
    FirstName = nm$first,
    MiddleName = nm$middle,
    LastName = nm$last,
    Suffix = "",
    Corresponding = if (isTRUE(a$corresponding)) "Yes" else "",
    HomePage = "",
    Group = "",
    ORCiD = a$orcid %||% "",
    stringsAsFactors = FALSE
  )
})
out <- do.call(rbind, rows)

# Use the template's own header, so column labels always match what bioRxiv wants.
header <- strsplit(readLines(TEMPLATE, n = 1, warn = FALSE), "\t")[[1]]
stopifnot(length(header) == ncol(out))
names(out) <- header

write.table(out, OUT, sep = "\t", quote = FALSE, row.names = FALSE, na = "")

message(sprintf("Wrote %s (%d authors)", OUT, nrow(out)))
missing_email <- out[[1]] == ""
if (any(missing_email)) {
  message(
    "Missing email for: ",
    paste(
      paste(out[[3]][missing_email], out[[5]][missing_email]),
      collapse = ", "
    )
  )
}
missing_orcid <- out[[10]] == ""
if (any(missing_orcid)) {
  message(
    "Missing ORCiD for: ",
    paste(
      paste(out[[3]][missing_orcid], out[[5]][missing_orcid]),
      collapse = ", "
    )
  )
}
