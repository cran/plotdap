% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/murSST.R
\docType{data}
\name{murSST}
\alias{murSST}
\title{murSST Data}
\format{
An object of class \code{griddap_nc} (inherits from \code{nc}, \code{data.frame}) with 0 rows and 2 columns.
}
\usage{
murSST
}
\description{
pre-Download of murSST in `add_griddap()` example so that example
can run within CRAN Time limits
}
\details{
obtained using the `rerddap` command
murSST <- griddap(
'jplMURSST41', latitude = c(22, 51), longitude = c(-140, -105),
time = c('last', 'last'), fields = 'analysed_sst'
)
}
\keyword{datasets}
