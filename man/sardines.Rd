% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sardines.R
\docType{data}
\name{sardines}
\alias{sardines}
\title{sardine Data}
\format{
An object of class \code{tabledap} (inherits from \code{data.frame}) with 56 rows and 5 columns.
}
\usage{
sardines
}
\description{
pre-Download of sardine data in `add_tabledap()` example so that example
can run within CRAN Time limits
}
\details{
obtained using the `rerddap` command
  sardines <- tabledap( 'FRDCPSTrawlLHHaulCatch',
     fields = c('latitude',  'longitude', 'time', '
         scientific_name', 'subsample_count'),
        'time>=2010-01-01', 'time<=2012-01-01',
        'scientific_name="Sardinops sagax"')
)
}
\keyword{datasets}
