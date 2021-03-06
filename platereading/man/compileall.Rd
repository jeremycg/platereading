% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/platereading.R
\name{compileall}
\alias{compileall}
\title{Filter a data frame by residuals, and combine passing values into means and sds}
\usage{
compileall(x, y = 0.01)
}
\arguments{
\item{x}{A dataframe containing the named output of the looper function}

\item{y}{A residual cutoff, defaults to 0.01}
}
\value{
A dataframe containing the data filtered and summarised
}
\description{
Filter a data frame by residuals, and combine passing values into means and sds
}

