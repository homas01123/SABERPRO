## code to prepare `a_w` dataset goes here
a_w <- read_csv(fs::path_package("SABER", "extdata", "a_w.csv"))
usethis::use_data(a_w, overwrite = TRUE)
