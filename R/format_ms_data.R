#' Tidy MS report
#'
#' format_ms_data() tidies Mascot search results provided by the Cambridge Centre for Proteomics (CCP) for further processing in R.
#' It takes a MS report saved as a txt file as input and returns a tibble summarising the results.
#'
#' @param filename a txt file derived from a mht MS data file. In Microsoft Word: Open the mht file, then "Save As..." > "File Format: Plain Text" > "Text encoding: Unicode UTF-8".
#'
#' @return A tibble summarising the Mascot search results.
#'
#' @export
#'
#' @importFrom dplyr filter
#' @importFrom dplyr group_by
#' @importFrom dplyr if_else
#' @importFrom dplyr lag
#' @importFrom dplyr lead
#' @importFrom dplyr mutate
#' @importFrom dplyr rename
#' @importFrom dplyr select
#' @importFrom dplyr ungroup
#' @importFrom magrittr %>%
#' @importFrom purrr map
#' @importFrom purrr map_chr
#' @importFrom purrr none
#' @importFrom purrr some
#' @importFrom readr read_delim
#' @importFrom stats setNames
#' @importFrom stringr str_detect
#' @importFrom stringr str_extract
#' @importFrom stringr str_remove
#' @importFrom tidyselect everything
#'
#' @examples
#'
format_ms_data <- function(filename) {

  ms_data <- read_delim(file = filename, delim = "\n") %>%
    setNames(nm = "raw_data") %>%
    mutate(raw_data = raw_data %>% str_remove("^\\s*") %>% str_remove("\\s*$")) %>%
    {.[-c(1:which(str_detect(.$raw_data, "^.*Error\\stolerant.*$"))), ]} %>%
    mutate(raw_data = if_else(str_detect(raw_data, "^\\s*$"), NA_character_, raw_data)) %>%
    filter(!is.na(raw_data)) %>%
    mutate(protein_hit_number = if_else(str_detect(raw_data, "^\\d*\\.\\s*$"), raw_data, NA_character_) %>% str_extract("\\d*"))

  cond <- some(ms_data$protein_hit_number, ~ is.na(.x))
  while(cond == TRUE) {
    ms_data <- ms_data %>% mutate(protein_hit_number = if_else(is.na(protein_hit_number), lag(protein_hit_number), protein_hit_number))
    cond <- some(ms_data$protein_hit_number, ~ is.na(.x))
  }

  ms_data <- ms_data %>%
    group_by(protein_hit_number) %>%
    nest %>%
    mutate(protein_id = map_chr(data, ~ .x[2,][[1]]),
           mass = protein_id %>% str_extract("Mass.*$") %>% str_remove("Mass:\\s") %>% str_remove("\\sS.*$") %>% as.double,
           score = protein_id %>% str_extract("Score.*$") %>% str_remove("Score:\\s") %>% str_remove("\\sM.*$") %>% as.double,
           matches = protein_id %>% str_extract("Matches.*$") %>% str_remove("Matches:\\s") %>% str_remove("\\sS.*$") %>% str_remove("\\s.*$"),
           sequences = protein_id %>% str_extract("Sequences.*$") %>% str_remove("Sequences:\\s") %>% str_remove("\\se.*$") %>% str_remove("\\s.*$"),
           emPAI = protein_id %>% str_extract("emPAI.*$") %>% str_remove("emPAI:\\s") %>% as.double,
           protein_id = protein_id %>% str_remove("^\\d::") %>% str_remove("\\s.*$"),
           protein_nm = map_chr(data, ~ .x[3,][[1]]),
           os = protein_nm %>% str_extract("OS=.*$") %>% str_remove("OS=") %>% str_remove("\\s[OX|GN|PE|SV].*$"),
           ox = protein_nm %>% str_extract("OX=.*$") %>% str_remove("OX=") %>% str_remove("\\s[GN|PE|SV].*$"),
           gn = protein_nm %>% str_extract("GN=.*$") %>% str_remove("GN=") %>% str_remove("\\s[PE|SV].*$"),
           pe = protein_nm %>% str_extract("PE=.*$") %>% str_remove("PE=") %>% str_remove("\\s.*$"),
           sv = protein_nm %>% str_extract("SV=.*$") %>% str_remove("SV=") %>% str_remove("\\s.*$"),
           protein_nm = protein_nm %>% str_remove("\\sOS.*$"),
           proteins_matching_same_peptide_set = map(data,
                                                    ~ if(none(.x$raw_data,
                                                                     ~ str_detect(.x, "Proteins matching the same set of peptides"))) {
                                                      NULL
                                                    } else {
                                                      .x[-c(1:which(str_detect(.x$raw_data, "Proteins matching the same set of peptides"))),] %>%
                                                        { if(some(.$raw_data,
                                                                         ~ str_detect(.x, "Proteins matching a subset of these peptides"))) {
                                                          .[-c(which(str_detect(.$raw_data, "Proteins matching a subset of these peptides")):nrow(.)),]
                                                        } else {
                                                          .
                                                        }
                                                        } %>%
                                                        rename(protein_id = raw_data) %>%
                                                        mutate(protein_nm = if_else(str_detect(protein_id, "^\\d::.*$"), NA_character_, protein_id),
                                                               protein_id = str_extract(protein_id, "^\\d::.*$") %>% str_remove("^\\d::"),
                                                               protein_nm = lead(protein_nm)) %>%
                                                        na.omit %>%
                                                        mutate(mass = protein_id %>% str_extract("Mass.*$") %>% str_remove("Mass:\\s") %>% str_remove("\\sS.*$") %>% as.double,
                                                               score = protein_id %>% str_extract("Score.*$") %>% str_remove("Score:\\s") %>% str_remove("\\sM.*$") %>% as.double,
                                                               matches = protein_id %>% str_extract("Matches.*$") %>% str_remove("Matches:\\s") %>% str_remove("\\sS.*$") %>% str_remove("\\s.*$"),
                                                               sequences = protein_id %>% str_extract("Sequences.*$") %>% str_remove("Sequences:\\s") %>% str_remove("\\se.*$") %>% str_remove("\\s.*$"),
                                                               emPAI = protein_id %>% str_extract("emPAI.*$") %>% str_remove("emPAI:\\s") %>% as.double,
                                                               protein_id = protein_id %>% str_remove("\\s.*$"),
                                                               os = protein_nm %>% str_extract("OS=.*$") %>% str_remove("OS=") %>% str_remove("\\s[OX|GN|PE|SV].*$"),
                                                               ox = protein_nm %>% str_extract("OX=.*$") %>% str_remove("OX=") %>% str_remove("\\s[GN|PE|SV].*$"),
                                                               gn = protein_nm %>% str_extract("GN=.*$") %>% str_remove("GN=") %>% str_remove("\\s[PE|SV].*$"),
                                                               pe = protein_nm %>% str_extract("PE=.*$") %>% str_remove("PE=") %>% str_remove("\\s.*$"),
                                                               sv = protein_nm %>% str_extract("SV=.*$") %>% str_remove("SV=") %>% str_remove("\\s.*$"),
                                                               protein_nm = protein_nm %>% str_remove("\\sOS.*$")) %>%
                                                        select(protein_id, protein_nm, everything())
                                                    }),
           proteins_matching_peptide_subset = map(data,
                                                  ~ if(none(.x$raw_data,
                                                                   ~ str_detect(.x, "Proteins matching a subset of these peptides"))) {
                                                    NULL
                                                  } else {
                                                    .x[-c(1:which(str_detect(.x$raw_data, "Proteins matching a subset of these peptides"))),] %>%
                                                      rename(protein_id = raw_data) %>%
                                                      mutate(protein_nm = if_else(str_detect(protein_id, "^\\d::.*$"), NA_character_, protein_id),
                                                             protein_id = str_extract(protein_id, "^\\d::.*$") %>% str_remove("^\\d::"),
                                                             protein_nm = lead(protein_nm)) %>%
                                                      na.omit %>%
                                                      mutate(mass = protein_id %>% str_extract("Mass.*$") %>% str_remove("Mass:\\s") %>% str_remove("\\sS.*$") %>% as.double,
                                                             score = protein_id %>% str_extract("Score.*$") %>% str_remove("Score:\\s") %>% str_remove("\\sM.*$") %>% as.double,
                                                             matches = protein_id %>% str_extract("Matches.*$") %>% str_remove("Matches:\\s") %>% str_remove("\\sS.*$") %>% str_remove("\\s.*$"),
                                                             sequences = protein_id %>% str_extract("Sequences.*$") %>% str_remove("Sequences:\\s") %>% str_remove("\\se.*$") %>% str_remove("\\s.*$"),
                                                             emPAI = protein_id %>% str_extract("emPAI.*$") %>% str_remove("emPAI:\\s") %>% as.double,
                                                             protein_id = protein_id %>% str_remove("\\s.*$"),
                                                             os = protein_nm %>% str_extract("OS=.*$") %>% str_remove("OS=") %>% str_remove("\\s[OX|GN|PE|SV].*$"),
                                                             ox = protein_nm %>% str_extract("OX=.*$") %>% str_remove("OX=") %>% str_remove("\\s[GN|PE|SV].*$"),
                                                             gn = protein_nm %>% str_extract("GN=.*$") %>% str_remove("GN=") %>% str_remove("\\s[PE|SV].*$"),
                                                             pe = protein_nm %>% str_extract("PE=.*$") %>% str_remove("PE=") %>% str_remove("\\s.*$"),
                                                             sv = protein_nm %>% str_extract("SV=.*$") %>% str_remove("SV=") %>% str_remove("\\s.*$"),
                                                             protein_nm = protein_nm %>% str_remove("\\sOS.*$")) %>%
                                                      select(protein_id, protein_nm, everything())
                                                  })) %>%
    select(protein_hit_number, protein_id, protein_nm,
           mass:emPAI, os:sv,
           proteins_matching_same_peptide_set,
           proteins_matching_peptide_subset,
           data) %>%
    ungroup()
}
