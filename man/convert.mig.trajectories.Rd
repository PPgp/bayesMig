% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predict_mig.R
\name{convert.mig.trajectories}
\alias{convert.mig.trajectories}
\title{Converting Trajectories of Migration Rates into ACSII Files}
\usage{
convert.mig.trajectories(
  sim.dir = NULL,
  n = 1000,
  output.dir = NULL,
  verbose = FALSE
)
}
\arguments{
\item{sim.dir}{Directory containing the prediction object. It should be the same as
the \code{output.dir} argument in \code{\link{mig.predict}}.}

\item{n}{Number of trajectories to be stored. It can be either a single number 
or the word \dQuote{all} in which case all available trajectories are converted.
If the number is smaller than the number of trajectories available 
in the prediction object, they are selected by equal spacing.}

\item{output.dir}{Directory into which the resulting files will be stored. 
If it is \code{NULL}, the same directory is used as for the prediction.}

\item{verbose}{Logical value. Switches log messages on and off.}
}
\value{
No return value.
}
\description{
Converts trajectories of the net migration rates stored 
    in a binary format into two CSV files.
}
\details{
The function creates two files. First, \dQuote{ascii_trajectories.csv}
    is a comma-separated table with the following columns: 
    \itemize{\item{\dQuote{LocID}: }{country code} 
        \item{\dQuote{Period}: }{prediction interval, e.g. 2015-2020} 
        \item{\dQuote{Year}: }{middle year of the prediction interval}
        \item{\dQuote{Trajectory}: }{identifier of the trajectory}
        \item{\dQuote{mig}: }{net migration rate}
    }
    The second file is called \dQuote{ascii_trajectories_wide.csv}, also 
    a comma-separated table and it contains the same information as above 
    but in a wide format. I.e. the data for one country are 
    ordered in columns, thus, there is one column per country. The country columns 
    are ordered alphabetically.
}
\note{
This function is automatically called from the \code{\link{mig.predict}} 
    function, therefore in standard cases it will not be needed to call it directly. 
    However, it can be useful for example, if different number of trajectories are to be converted, 
    without having to re-run the prediction.
}
\seealso{
\code{\link[bayesTFR]{convert.tfr.trajectories}}, 
  \code{\link{mig.write.projection.summary}}, \code{\link{get.mig.trajectories}}
}
