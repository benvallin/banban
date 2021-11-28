library(tidyverse)
library(xml2)

import_cell_profiler_data <- function(RelateObjects.data.files.list,
                                      IdentifyPrimaryObjects.data.files.list,
                                      prefix,
                                      main.object.nm = "neuron",
                                      main.object.structure.nm = "cell_body",
                                      metadata = FALSE) {

  message("\nProcessing CellProfiler data...\n")

  # Define the main object name and assign it to a variable, by default it's neuron -> "neuron_id"
  main.object.var.nm <- str_c(main.object.nm, "_id")

  # create a dataframe RelateObjects.data (2,6): first col contains the name of the object, second col contains the relate object data from cell profiler, usually 6 rows
  # 1 - create RelateObjects.data which contains the name of the relate object csv files (1,6), col called "RelateObjects.nm" contains the names include in the list.txt created with bash
  # 2 - map along the first col to create a second one (2,6): RelateObjects.val which contains the actual data + renaming FileName -> image_id and Imagenumber -> image_number
  # 3 - mutate the first col to get a meaningful name without RelateObjects_ and csv
  RelateObjects.data <- read_delim(str_c(prefix, RelateObjects.data.files.list),
                                   col_names = "RelateObjects.nm",
                                   delim = "\n",
                                   show_col_types = FALSE) %>%
    mutate(RelateObjects.val = map(RelateObjects.nm,
                                   ~ read_csv(str_c(prefix, .x), show_col_types = FALSE, progress = FALSE) %>%
                                     mutate(image_id = str_remove(FileName_czi, ".czi"))),
           RelateObjects.nm = str_remove(RelateObjects.nm, "RelateObjects_") %>% str_remove(., ".csv"))

  # create a dataframe main.object.data (2,1) and create a col named main.object.var.nm containing "img_x_neuron_x" that we will need to join later
  # 1 - from RelateObjects.data, keep only the row which match value of the main.object.structure.nm, by default "cell_body", table is now (2,1)
  # 2 - unest the content of the second col, RelateObjects.val
  # 3 - create a col main.object.var.nm that concatenate image number and object number to give for example "img_2_neuron_3" -- not familliar with quasiquotation AMOUR
  # 4 - tidy the table by putting meaningful col first + the rest and removing the useless one -- col image_number is select then deleted AMOUR
  # 5 - keep the col we wanted to create main.object.var.nm that contains "img_x_neuron_x" + nest all the rest into one = "cell_body" + "_data"
  main.object.data <- RelateObjects.data %>%
    filter(RelateObjects.nm == main.object.structure.nm) %>%
    unnest(cols = RelateObjects.val) %>%
    mutate(!!main.object.var.nm := str_c(image_id, "_", main.object.nm, "_", ObjectNumber)) %>%
    select(image_id, all_of(main.object.var.nm), everything()) %>%
    nest(!!str_c(main.object.structure.nm, "_data") := -all_of(main.object.var.nm))

  # create a dataframe IdentifyPrimaryObjects.data
  # 1 - create IdentifyPrimaryObjects.data which contains the name of the Identify primary object csv files,
  #     first col called "IdentifyPrimaryObjects.nm" contains the names include in the list.txt created with bash
  # 2 - map along the first col to create a second one: IdentifyPrimaryObjects.val which contains the actual data
  # 3 - mutate the first col to get a meaningful name without "IdentifyPrimaryObjects_" and csv
  # 4 - map logical along the second col to create has.Parent_RelateObjects_main.object.structure.nm,
  #     map check in IdentifyPrimaryObjects.val if there is a col named "Parent_RelateObjects_" + "cell_body" (by default) and return True or False
  # 5 - filter the data frame by removing the row in which the value in the col "has.Parent_RelateObjects_main.object.structure.nm" is false
  #     discard all the object that do not have parent, meaning they do not belong to a cell body
  # 6 - mutate second col: go through both IdentifyPrimaryObjects.val and IdentifyPrimaryObjects.nm
  #     along IdentifyPrimaryObjects.val, change ImageNumber in image_number + ObjectNumber in object_number, create main.object_number = AMOUR
  #     filter rows that have a main object number different from 0
  #     create main.object.var.nm containing "img_x_neuron_x",
  #     create temp_object_id by going through IdentifyPrimaryObjects.nm and adding object number, create object type by going through IdentifyPrimaryObjects.nm
  IdentifyPrimaryObjects.data <- read_delim(str_c(prefix, IdentifyPrimaryObjects.data.files.list),
                                            col_names = "IdentifyPrimaryObjects.nm",
                                            delim = "\n",
                                            show_col_types = FALSE) %>%
    mutate(IdentifyPrimaryObjects.val = map(IdentifyPrimaryObjects.nm, ~ read_csv(str_c(prefix, .x), show_col_types = FALSE, progress = FALSE)),
           IdentifyPrimaryObjects.nm = str_remove(IdentifyPrimaryObjects.nm, "IdentifyPrimaryObjects_") %>% str_remove(".csv"),
           has.Parent_RelateObjects_main.object.structure.nm = map_lgl(IdentifyPrimaryObjects.val,
                                                                       ~ if_else(some(str_detect(colnames(.x), str_c("Parent_RelateObjects_", main.object.structure.nm)), isTRUE),
                                                                                 TRUE,
                                                                                 FALSE))) %>%
    filter(has.Parent_RelateObjects_main.object.structure.nm == TRUE) %>%
    select(-has.Parent_RelateObjects_main.object.structure.nm) %>%
    mutate(IdentifyPrimaryObjects.val = map2(IdentifyPrimaryObjects.val, IdentifyPrimaryObjects.nm,
                                             ~ .x %>%
                                               mutate(image_id = str_remove(FileName_czi, ".czi"),
                                                      main.object_number = .[[str_c("Parent_RelateObjects_", main.object.structure.nm)]]) %>%
                                               filter(main.object_number != 0) %>%
                                               mutate(!!main.object.var.nm := str_c(image_id, "_", main.object.nm, "_", main.object_number),
                                                      temp_object_id = str_c(image_id, "_", .y, "_", ObjectNumber),
                                                      object_type = .y) %>%
                                               nest(IdentifyPrimaryObjects.val_data = c(ObjectNumber, Number_Object_Number))),
           IdentifyPrimaryObjects.val.empty = map_lgl(.x = IdentifyPrimaryObjects.val,
                                                      .f = ~ifelse(nrow(.x) == 0L, TRUE, FALSE))) %>%
    filter(IdentifyPrimaryObjects.val.empty == FALSE) %>%
    select(-IdentifyPrimaryObjects.val.empty)

  internal.object.data <- RelateObjects.data %>%
    filter(RelateObjects.nm %in% IdentifyPrimaryObjects.data$IdentifyPrimaryObjects.nm) %>%
    mutate(RelateObjects.val = map2(RelateObjects.val, RelateObjects.nm,
                                    ~ .x %>%
                                      mutate(temp_object_id = str_c(image_id, "_", .y, "_", .[[str_c("Parent_IdentifyPrimaryObjects_", .y)]]),
                                             object_type = .y) %>%
                                      nest(RelateObjects.val_data = c(ObjectNumber, Number_Object_Number))))

  all.data <- map2(IdentifyPrimaryObjects.data$IdentifyPrimaryObjects.val, internal.object.data$RelateObjects.val,
                   ~ full_join(.x, .y,
                               by = c("ImageNumber", "FileName_czi", "PathName_czi", "Location_Center_X", "Location_Center_Y", "Location_Center_Z", "image_id", "temp_object_id", "object_type")) %>%
                     mutate(object_id = str_c(.[[!!main.object.var.nm]], "_", temp_object_id)) %>%
                     select(object_id, everything())
                   ) %>%
    map(~ nest(.data = .x, data = -c(all_of(main.object.var.nm), object_type))) %>%
    bind_rows() %>%
    pivot_wider(names_from = "object_type", values_from = "data") %>%
    rename_if(.predicate = function(x) is.list(x), .funs = function(x) str_c(x, "_data")) %>%
    full_join(main.object.data, , by = str_c(main.object.nm, "_id")) %>%
    mutate(image_id = map_chr(.[[str_c(main.object.structure.nm, "_data")]], ~ .x$image_id),
           file_id = map_chr(.[[str_c(main.object.structure.nm, "_data")]], ~ .x$FileName_czi)) %>%
    select(file_id, image_id, all_of(main.object.var.nm), str_c(main.object.structure.nm, "_data"), everything())

  if(metadata == TRUE) {

    message("Processing metadata...\n")

    all.data <- full_join(all.data,
                          all.data[,"image_id"] %>%
                            unique %>%
                            mutate(metadata = map(image_id,
                                                  ~ suppressMessages(read_xml(str_c(prefix, "metadata/", .x,".tif_metadata.xml"),
                                                                              options = c("NOBLANKS")) %>%
                                                                       as_list() %>%
                                                                       {.$ImageMetadata$Scaling$Items[1:2]} %>%
                                                                       map(~ as.double(.x$Value[[1]])) %>%
                                                                       bind_cols %>%
                                                                       setNames(nm = c("scaling_x", "scaling_y"))))) %>%
                            unnest(cols = metadata),
                          by = "image_id")
  }

  message("Done!\n")

  all.data

}


test <- import_cell_profiler_data(RelateObjects.data.files.list = "list.relate.objects.files.txt",
                                  IdentifyPrimaryObjects.data.files.list = "list.primary.objects.files.txt",
                                  prefix = "~/Desktop/cp files/doesn't work/",
                                  main.object.nm = "neuron",
                                  main.object.structure.nm = "cell_body",
                                  metadata = FALSE)


test2 <- banban::import_cell_profiler_data(RelateObjects.data.files.list = "list.relate.objects.files.txt",
                                           IdentifyPrimaryObjects.data.files.list = "list.primary.objects.files.txt",
                                           prefix = "~/Desktop/cp files/doesn't work/",
                                           main.object.nm = "neuron",
                                           main.object.structure.nm = "cell_body",
                                           metadata = FALSE)

x <- test$lc3_data[[426]]$AreaShape_Area

y <- test2$lc3_data[[426]]$AreaShape_Area
