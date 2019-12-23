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
library(tmsClassifier)  
  
n.subjects <- 10  
TMS <- tms.coef[sample(nrow(tms.coef), n.subjects), ]  
predicted <- tmsClassify(TMS)
```

In this example, we used estimated TMS regression coefficients from default data (object `tms.coef`). The user may generate a custom TMS coefficients matrix using `tmsClassifier` functions (see package documentation for details).

## Input TMS data

Package tmsClassifier uses a default TMS data format. Subject classification is based on three blocks of attributes: (i) TMS measures, (ii) diagnosis field (optional), (iii) covariates (optional).  
The default dataset includes 4 TMS measures:
- SICI-ICF (SICI: short-interval intracortical inhibition; ICF: intracortical facilitation) as the first 7 columns, taken at times (interstimulus intervals, ISI): 1, 2, 3, 5 ms (SICI) and 7, 10, 15 ms (ICF);
- SAI (short-latency afferent inhibition). By default, they should be the 4 columns following SICI-ICF, taken at time steps (ISI): -4, 0, 4, 8 ms;
- LICI (long-interval intracortical inhibition). By default, they should be the 3 columns following SAI, taken at time
steps (ISI): 50, 100, 150 ms.

```{r, eval = FALSE}
  SICI_1  SICI_2 SICI_3 SICI_5     ICF_7 ICF_10   ICF_15    SAI_m4   SAI_0
1  0.405 0.63600  0.797  1.510 1.9640000  2.057 1.390000 0.7800000 0.92500
2  0.061 0.05800  0.100  1.159 1.1810000  1.745 1.694000 0.8500000 0.87400
3  0.570 0.72500  0.699  1.085 1.8420000  1.875 2.418000 1.1320000 0.47900
4  0.000 0.00000  0.444  0.556 0.8847778  1.111 1.262250 0.9270000 0.54225
5  0.080 0.19425  0.150  0.170 1.1300000  1.300 1.301875 0.8518571 0.42175
6  1.031 0.95900  0.791  1.082 1.1480000  0.843 0.780000 0.9420000 0.38600
7  0.406 0.80800  1.492  1.649 1.5080000  1.001 2.126000 0.9640000 1.07400
8  1.142 1.87200  1.872  1.121 1.7080000  1.738 1.766000 0.9255000 1.06100
9  0.097 0.60300  0.851  0.837 0.9810000  1.002 0.920000 0.8690000 0.60000
10 0.380 0.78100  1.081  1.404 1.5190000  1.311 1.387000 0.9630000 1.28600
   SAI_4 SAI_8   LICI_50  LICI_100  LICI_150 diagnosis center sex age
1  0.772 0.557 0.1590000 0.4190000 0.4250000        AD     C1   M  76
2  0.779 0.899 0.2730000 0.3872857 0.4587143        AD     C1   F  68
3  0.283 0.592 0.3397143 0.6042857 0.5590000        HC     C3   F  63
4  0.673 0.577 0.2573333 0.1993333 0.5286667        HC     C2   M  76
5  0.450 0.450 0.2380000 0.1784000 0.3358000       DLB     C4   F  77
6  0.563 1.070 0.0300000 0.1560000 0.6860000       FTD     C1   F  63
7  0.663 0.888 0.4985000 0.7678750 0.5547500        AD     C1   M  57
8  1.053 1.038 0.2921429 0.4844286 0.4631429        AD     C2   F  68
9  0.489 0.858 0.5810000 0.6040000 0.4580000       FTD     C1   M  52
10 0.860 0.890 0.5860000 0.5800000 0.5900000       DLB     C1   M  78
```


