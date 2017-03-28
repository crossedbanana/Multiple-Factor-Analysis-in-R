#workflow

library(devtools)
devtools::document()
# checking documentation
devtools::check_man()
# running tests
devtools::test()
devtools::build_vignettes()
# building tarball (e.g. oski_0.1.tar.gz)
devtools::build()
# checking install
devtools::install()

