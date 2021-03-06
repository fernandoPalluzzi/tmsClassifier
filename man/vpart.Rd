% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TMS_Classifier.R
\name{vpart}
\alias{vpart}
\title{Data partitioning utility}
\usage{
vpart(data, p, status = NULL, shuffle = FALSE)
}
\arguments{
\item{data}{An input data.frame.}

\item{p}{Proportion of input rows to be extracted from the input
dataset.}

\item{status}{Numeric value indicating the column for subject diagnosis
(default status = NULL). Allowed values for this column are "HC", "FTD",
"AD", and "DLB". If status is not NULL, random partitioning will preserve
the status column proportions.}

\item{shuffle}{A logical value. If TRUE, the input rows are randomly
shuffled before data partitioning.}
}
\value{
A list of 2 data.frames:
\enumerate{
\item "training.set", the portion of the input data defined by p;
\item "validation.set", the portion of the input data defined by 1-p.
}
}
\description{
Extract a random partition from an input dataset.
}
