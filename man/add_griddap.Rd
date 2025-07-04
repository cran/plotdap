% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/add_griddap.R
\name{add_griddap}
\alias{add_griddap}
\title{Add rerddap::griddap() data to a plotdap map}
\usage{
add_griddap(
  plot,
  grid,
  var,
  fill = "viridis",
  maxpixels = 10000,
  time = mean,
  animate = FALSE,
  cumulative = FALSE,
  ...
)
}
\arguments{
\item{plot}{a \link{plotdap} object.}

\item{grid}{a \link[rerddap]{griddap} object.}

\item{var}{a formula defining a variable, or function of variables to visualize.}

\item{fill}{either a character string of length 1 matching a name in the
package \code{cmocean} or a vector of color codes.
This defines the colorscale used to encode values
of \code{var}.}

\item{maxpixels}{integer > 0. Maximum number of cells to use for the plot.
If maxpixels < ncell(x), sampleRegular is used before plotting.
If gridded=TRUE maxpixels may be ignored to get a larger sample}

\item{time}{how to resolve multiple time frames. Choose one of the following:
\itemize{
  \item A function to apply to each observation at a particular location
  (\link{mean} is the default).
  \item A character string (of length 1) matching a time value.
}}

\item{animate}{whether to animate over the \code{time} variable (if it exists).
Currently only implemented for \code{method='ggplot2'} and requires the
gganimate package.}

\item{cumulative}{- if animation should be cumulative -default FALSE}

\item{...}{arguments passed along to \code{geom_sf()}
(if \code{method='ggplot2'}, otherwise ignored).}
}
\value{
A plotdap object
}
\description{
\code{add_griddap} adds the data from an 'rerddap::griddap() call to
a 'plotdap' map
}
\examples{

# base plotting tends to be faster,
# but is less extensible plotdap("base")

# actual datasets in data folder to meet execution timings


 # murSST <- rerddap::griddap(
 #  ' jplMURSST41', latitude = c(35, 40), longitude = c(-125, -120.5),
 #   time = c('last', 'last'), fields = 'analysed_sst'
 # )

 # QMwind <- rerddap::griddap(
 #  'erdQMwindmday', time = c('2016-11-16', '2017-01-16'),
 #  latitude = c(30, 50), longitude = c(210, 240),
 #  fields = 'x_wind'
 #  )

p <- plotdap(crs = "+proj=robin")
p <- add_griddap(p, murSST, ~analysed_sst)

 # p <- plotdap(mapTitle = "Average wind over time")
 # p <- add_griddap(p, QMwind, ~x_wind)

# p <- plotdap("base", crs = "+proj=robin")
# p <- add_griddap(p, murSST, ~analysed_sst)

# layer tables on top of grids
require(magrittr)
p <- plotdap("base") \%>\%
  add_griddap(murSST, ~analysed_sst) \%>\%
  add_tabledap(sardines, ~subsample_count)

# multiple time periods
p <- plotdap("base", mapTitle = "Average wind over time")
p <- add_griddap(p, QMwind, ~x_wind)

}
