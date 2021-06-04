#' Import CellProfiler data
#'
#' blabla
#'
#' @param RelateObjects.data.files.list blabla
#' @param IdentifyPrimaryObjects.data.files.list blabla
#' @param prefix blabla
#' @param main.object.nm blabla
#' @param main.object.structure.nm blabla
#' @param metadata blabla
#'
#' @return blabla
#'
#' @export
#'
#' @importFrom dplyr bind_rows
#' @importFrom dplyr filter
#' @importFrom dplyr full_join
#' @importFrom dplyr if_else
#' @importFrom dplyr mutate
#' @importFrom dplyr rename
#' @importFrom dplyr rename_if
#' @importFrom dplyr select
#' @importFrom magrittr %>%
#' @importFrom purrr map
#' @importFrom purrr map_chr
#' @importFrom purrr map_lgl
#' @importFrom purrr map2
#' @importFrom purrr some
#' @importFrom readr read_csv
#' @importFrom readr read_delim
#' @importFrom rlang !!
#' @importFrom rlang :=
#' @importFrom rlang as_list
#' @importFrom stats setNames
#' @importFrom stringr str_c
#' @importFrom stringr str_detect
#' @importFrom stringr str_remove
#' @importFrom tidyr nest
#' @importFrom tidyr pivot_wider
#' @importFrom tidyr unnest
#' @importFrom tidyselect everything
#' @importFrom xml2 read_xml
#'
#' @examples
#'
import_cell_profiler_data <- function(RelateObjects.data.files.list,
                                      IdentifyPrimaryObjects.data.files.list,
                                      prefix,
                                      main.object.nm = "neuron",
                                      main.object.structure.nm = "cell_body",
                                      metadata = FALSE) {

  main.object.var.nm <- str_c(main.object.nm, "_id")

  RelateObjects.data <- read_delim(str_c(prefix, RelateObjects.data.files.list), col_names = "RelateObjects.nm", delim = "\n") %>%
    mutate(RelateObjects.val = map(RelateObjects.nm,
                                   ~ read_csv(str_c(prefix, .x)) %>%
                                     mutate(image_id = str_remove(FileName_czi, ".czi")) %>%
                                     rename(image_number = ImageNumber))) %>%
    mutate(RelateObjects.nm = str_remove(RelateObjects.nm, "RelateObjects_") %>% str_remove(., ".csv"))

  main.object.data <- RelateObjects.data %>%
    filter(RelateObjects.nm == main.object.structure.nm) %>%
    unnest(cols = RelateObjects.val) %>%
    mutate(!!main.object.var.nm := str_c("img_", image_number, "_", main.object.nm, "_", ObjectNumber)) %>%
    select(image_id, image_number, main.object.var.nm, everything(), -c(image_number, RelateObjects.nm, FileName_czi, PathName_czi)) %>%
    nest(!!str_c(main.object.structure.nm, "_data") := -main.object.var.nm)

  IdentifyPrimaryObjects.data <- read_delim(str_c(prefix, IdentifyPrimaryObjects.data.files.list), col_names = "IdentifyPrimaryObjects.nm", delim = "\n") %>%
    mutate(IdentifyPrimaryObjects.val = map(IdentifyPrimaryObjects.nm, ~ read_csv(str_c(prefix, .x))),
           IdentifyPrimaryObjects.nm = str_remove(IdentifyPrimaryObjects.nm, "IdentifyPrimaryObjects_") %>% str_remove(".csv"),
           has.Parent_RelateObjects_main.object.structure.nm = map_lgl(IdentifyPrimaryObjects.val,
                                                                       ~ if_else(some(str_detect(colnames(.x), str_c("Parent_RelateObjects_", main.object.structure.nm)), isTRUE),
                                                                                 TRUE,
                                                                                 FALSE))) %>%
    filter(has.Parent_RelateObjects_main.object.structure.nm == TRUE) %>%
    mutate(IdentifyPrimaryObjects.val = map2(IdentifyPrimaryObjects.val, IdentifyPrimaryObjects.nm,
                                             ~ .x %>% select(image_number = ImageNumber,
                                                             main.object_number = !!str_c("Parent_RelateObjects_", main.object.structure.nm),
                                                             object_number = ObjectNumber) %>%
                                               filter(main.object_number != 0) %>%
                                               mutate(!!main.object.var.nm := str_c("img_", image_number, "_", main.object.nm, "_", main.object_number),
                                                      temp_object_id = str_c(.y, "_", object_number),
                                                      object_type = .y) %>%
                                               select(-c(main.object_number, object_number))
    )) %>%
    select(-has.Parent_RelateObjects_main.object.structure.nm)

  internal.object.data <- RelateObjects.data %>%
    filter(RelateObjects.nm %in% IdentifyPrimaryObjects.data$IdentifyPrimaryObjects.nm) %>%
    mutate(RelateObjects.val = map2(RelateObjects.val, RelateObjects.nm,
                                    ~ rename(.x, object_number = str_c("Parent_IdentifyPrimaryObjects_", .y)) %>%
                                      mutate(temp_object_id = str_c(.y, "_", object_number),
                                             object_type = .y) %>%
                                      select(image_id, image_number, temp_object_id, everything(), -object_number)))

  all.data <- map2(IdentifyPrimaryObjects.data$IdentifyPrimaryObjects.val, internal.object.data$RelateObjects.val,
                   ~ full_join(.x, .y) %>%
                     mutate(object_id = str_c(.[[!!main.object.var.nm]], "_", temp_object_id)) %>%
                     select(object_id, everything(), -c(image_number, temp_object_id, image_id, FileName_czi, PathName_czi))) %>%
    map(~ nest(.data = .x, data = -c(main.object.var.nm, object_type))) %>%
    bind_rows() %>%
    pivot_wider(names_from = "object_type", values_from = "data") %>%
    rename_if(.predicate = function(x) is.list(x), .funs = function(x) str_c(x, "_data")) %>%
    full_join(main.object.data) %>%
    mutate(image_id = map_chr(.[[str_c(main.object.structure.nm, "_data")]], ~ .x$image_id),
           !!str_c(main.object.structure.nm, "_data") := map(.[[str_c(main.object.structure.nm, "_data")]], ~ select(.x, -image_id))) %>%
    select(image_id, main.object.var.nm, str_c(main.object.structure.nm, "_data"), everything()) %>%
    {if(metadata == TRUE) {
      full_join(., .[,"image_id"] %>%
                  unique %>%
                  mutate(metadata = map(image_id, ~ read_xml(str_c(prefix, "metadata/", .,".tif_metadata.xml")) %>%
                                          as_list() %>%
                                          {.$ImageMetadata$Scaling$Items[1:2]} %>%
                                          map(~ as.double(.x$Value[[1]])) %>%
                                          bind_cols %>%
                                          setNames(nm = c("scaling_x", "scaling_y")))) %>%
                  unnest(cols = metadata))
    } else {
      .
    }}
}
