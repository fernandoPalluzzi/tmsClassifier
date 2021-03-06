% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TMS_Classifier.R
\name{tmsClassify}
\alias{tmsClassify}
\title{TMS Random Forest Classifier (RFC)}
\usage{
tmsClassify(data, cv = FALSE, status = 10, p = 1, k = 5)
}
\arguments{
\item{data}{A data.frame containing subjects as rows and TMS parameters
as columns. Optionally, a column specifying patient diagnosis can be
added to test RFC performances using cross-validation.}

\item{cv}{A logical value. If cv is TRUE, the cross-validation mode is
enabled, and a status (i.e., diagnosis) column is needed (see below).
By default, cv = FALSE and the diagnosis prediction mode is enabled.}

\item{status}{Numeric value indicating the column for subject diagnosis
(default status = 10). This value is ignored if cv is set to FALSE.}

\item{p}{Proportion of the input dataset to be used for the
cross-validation. For instance, if p = 0.75, a randomly chosen 25% of
the input data.frame will be isolated and returned as validation set,
while the remaining 75% is used for the cross-validation. Every partition
is balanced so that the original class proportion is preserved. By default,
p = 1 (i.e., the entire dataset is used for cross-validation) and no
validation set is generated.}

\item{k}{Number of cross-validation iterations (default k = 5).}
}
\value{
If cv = FALSE, the vector of diagnosis prediction for each
subject is returned. If cv = TRUE, an object of class cross.validation
is generated, including cross-validation results. If p < 1, a list of
three objects is returned:
\enumerate{
\item "CV", cross-validation performances;
\item "trainingSet", dataset partition used for the cross-vaidation;
\item "validationSet", dataset partition isolated before cross-validation.
}
}
\description{
Classify subjects with unknown diagnosis in one class among
Healthy subject (HC), Frontotemporal Dementia (FTD), Alzheimer's Disease
(AD), or Dementia with Lewy Bodies (DLB), using estimated TMS regression
parameters. This function provides a default RFC, requiring an input
data.frame of 9 columns (4 for SICI-ICF, 3 for SAI, and 2 for LICI; see
\code{\link[tmsClassifier]{tmsRegression}} for details). See function
\code{\link[tmsClassifier]{buildPredictor}} to build a custom RFC.
Alternatively, if a diagnosis column is present in the input data.frame,
it can be used to test RFC performances through cross-validation.
}
\references{
Slawski M, Daumer M, Boulesteix AL (2008). CMA - a comprehensive 
Bioconductor package for supervised classification with high 
dimensional data. BMC Bioinformatics, 9:439. 
https://doi.org/10.1186/1471-2105-9-439
}
\seealso{
\code{\link[CMA]{GenerateLearningsets}} and
\code{\link[CMA]{classification}}.
}
