# tmsClassifier
Random Forest Classifiers for Transcranial Magnetic Stimulation Data

The `tmsClassifier` package predicts subject diagnosis (Frontotemporal Dementia, Alzheimer's Disease, Dementia with Lewy Bodies, or Healthy) based on Transcranial Magnetic Stimulation (TMS) data, through the integration of linear modeling and Random Forest Classifiers.

## Installation

The development version of `tmsClassifier` can be installed from GitHub:

```{r, eval = FALSE}
devtools::install_github("fernandoPalluzzi/tmsClassifier")
```

## Unknown subject classification

With the term of "unknown subject" we define a patient with no ascertained diagnosis. Package `tmsClassifier` provides a preset TMS coefficient matrix, that can be used to predict unknown subject diagnosis.
Default TMS signatures include: Alzheimer's Disease (AD), Frontotemporal Dementia (FTD), Dementia with Lewy Bodies (DLB), and healthy controls (HC).
The prediction can be simply done with the `tmsClassify` command:

```{r, eval = FALSE}
n.subjects <- 10  
TMS <- tms.coef[sample(nrow(tms.coef), n.subjects), ]  
predicted <- tmsClassify(TMS)
```

In this example, we used estimated TMS regression coefficients from default data (object `tms.coef`). The user may generate a custom TMS coefficients matrix using `tmsClassifier` functions (see package documentation for details).
