## code to prepare `a_star_chl` dataset goes here
library(readr)
a_star_chl <- read_delim(
  fs::path_package("SABER", "extdata", "astarchl.txt"),
  delim = "		",
  skip = 9,
  n_max = 91,
  col_names = c("wavelength", "a_star_chl")
  )

usethis::use_data(a_star_chl, overwrite = TRUE)
