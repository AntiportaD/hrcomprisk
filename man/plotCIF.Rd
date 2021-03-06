% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotCIF.R
\name{plotCIF}
\alias{plotCIF}
\title{Plot Incidence and Ratio of sHR/csHR}
\usage{
plotCIF(cifobj, maxtime = Inf, ci = NULL, eoi = -1)
}
\arguments{
\item{cifobj}{A dataframe containing the Cumulative Incidence of each competing event by exposure group.}

\item{maxtime}{Largest time to display on the x-axis of all output plots.}

\item{ci}{A dataframe containing the 95\% CI for each ratio of csHR/sHR.}

\item{eoi}{Event number for the event of interest, useful when more than two events exist.}
}
\value{
A large list containing 2 figures:
\describe{
  \item{$plot1 }{  Plot the cumulative incidence of the composite event and of each event by exposure group.}
  \item{$plot2 }{  Plot the ratio of Sub-hazard ratio and Cause-specific hazard ratio for each event i (Ri).}
}
}
\description{
Plot Cumulative Incidence and Ratio of sHR/csHR.
}
\examples{
#data from the package
data <- hrcomprisk::dat_ckid
#Estimate the Cumulative Incidence Functions and Ratios of sHR and csHR
mydat.CIF<-CRCumInc(df=data, time=exit, event=event, exposed=b1nb0, print.attr=FALSE)
#Plot the CIs and Ratios estimated
plots<-plotCIF(cifobj=mydat.CIF, maxtime = 20, eoi = 1)
}
