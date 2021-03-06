% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nlme_fit_prep.R
\name{removeFailedDrugs}
\alias{removeFailedDrugs}
\title{Removes well positions where the drug is now considered unsuitable for
 screening from GDSC raw data, i.e., TAG = 'FAIL'}
\usage{
removeFailedDrugs(myDat)
}
\arguments{
\item{myDat}{a GDSC raw data data frame.}
}
\description{
\code{removeFailedDrugs} removes rows from GDSC raw data where the
 \code{TAG == 'FAIL'}.
}
\examples{
data("gdsc_example")
gdsc_example <- removeFailedDrugs(gdsc_example)
gdsc_example <- removeMissingDrugs(gdsc_example)
gdsc_example <- normalizeData(gdsc_example)
gdsc_example <- setConcsForNlme(gdsc_example)
nlme_data <- prepNlmeData(gdsc_example, "COSMIC_ID")
 
}
\seealso{
\code{\link{removeMissingDrugs}}, \code{\link{normalizeData}},
 \code{\link{setConcsForNlme}}, \code{\link{prepNlmeData}}
}
