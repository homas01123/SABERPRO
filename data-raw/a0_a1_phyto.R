## code to prepare `a0_a1_phyto` dataset goes here
a0_a1_phyto <- read_csv(fs::path_package("SABER", "extdata", "a0_a1_phyto.csv"))
usethis::use_data(a0_a1_phyto, overwrite = TRUE)
