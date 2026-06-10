## code to prepare `nomad` dataset goes here
nomad <- read_csv(fs::path_package("SABER", "extdata", "nomad_simplified.csv"))
usethis::use_data(nomad, overwrite = TRUE)
