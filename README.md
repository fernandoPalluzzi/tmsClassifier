# tmsClassifier
Random Forest Classifiers for Transcranial Magnetic Stimulation Data

The **tmsClassifier** package predicts subject diagnosis (Frontotemporal Dementia, Alzheimer's Disease, Dementia with Lewy Bodies, or Healthy) based on Transcranial Magnetic Stimulation (TMS) data, through the integration of linear modeling and Random Forest Classifiers.

## Installation

The development version of **tmsClassifier** can be installed from GitHub:

```{r, eval = FALSE}
devtools::install_github("fernandoPalluzzi/tmsClassifier")
```

## Unknown subject classification

With the term of "unknown subject" we define a patient with no ascertained diagnosis. Package **tmsClassifier** provides a default TMS coefficient matrix, that can be used to predict unknown subject diagnosis.
Default TMS signatures include: Alzheimer's Disease (AD), Frontotemporal Dementia (FTD), Dementia with Lewy Bodies (DLB), and healthy controls (HC).
The prediction can be simply done with the `tmsClassify` command:

```{r, eval = FALSE}
library(tmsClassifier)

# Sample 10 random subjects from default TMS data
n.subjects <- 10
TMS <- tms.coef[sample(nrow(tms.coef), n.subjects), ]

# Subjects' diagnosis prediction
predicted <- tmsClassify(TMS)
```

The input `TMS` object must be a data.frame containing TMS regression coefficients as columns and subjects as rows. In this example, we used estimated TMS regression coefficients from default data (object `tms.coef`). The user may generate a custom TMS coefficients matrix using **tmsClassifier** functions (see below and package documentation for details).

## Input TMS data and imputation

The TMS regression coefficient matrix can be directly generated from TMS data and, optionally, subject covariates (e.g., age and sex). Package **tmsClassifier** uses a default TMS data format, including three blocks of attributes (columns): (i) TMS measures, (ii) diagnosis field (optional), and (iii) covariates (optional).  
The default dataset includes 4 TMS measures:
- SICI-ICF (SICI: short-interval intracortical inhibition; ICF: intracortical facilitation) as the first 7 columns, taken at times (interstimulus intervals, ISI): 1, 2, 3, 5 ms (SICI) and 7, 10, 15 ms (ICF);
- SAI (short-latency afferent inhibition). By default, they should be the 4 columns following SICI-ICF, taken at time steps (ISI): -4, 0, 4, 8 ms;
- LICI (long-interval intracortical inhibition). By default, they should be the 3 columns following SAI, taken at time
steps (ISI): 50, 100, 150 ms.  

Optionally, the user may provide a column specifying subjects' diagnosis (e.g., the `diagnosis` field in the example below). This field is not required for diagnosis prediction, but only to test RFC performances (see **tmsClassifier** documentation). In the TMS sample data below, the last three columns are common covariates: `center` (i.e., where data were generated), `sex` (i.e., subject's sex), and `age` (i.e., subject's age at which TMS was done). If required by the user, TMS regression parameters can be adjusted to remove covariates influence.

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

If TMS data contain missing values (`NA`), they can be inputed with the `tmsImpute` command:

```{r, eval = FALSE}
tms <- tmsImpute(tms.withNA)
```

yielding the input dataset with imputed missing values. 

##  TMS regression parameters

If TMS data do not contain missing values, it can be used to generate TMS regression coefficients. This matrix is used as input for subject diagnosis prediction and other **tmsClassifier** functions. For each subject, TMS indicators (SICI-ICF, SAI, LICI) are modeled as polynomial functions of time (t), in the form y ~ poly(t). This function estimates two parameters for SICI (SICI = bs0 + bs\*t), two parameters for ICF (ICF = bi0 + bi\*t), three parameters for SAI (SAI = b0 + b1\*t + b2\*t^2), and two parameters for LICI (LICI = a0 + a1\*t). Regression coefficients can be generated using the command:

```{r, eval = FALSE}
# No covariates adjustment
B <- tmsRegression(tms)

# Adjustment by center (column 16), sex (column 17), and age (column 18)
B <- tmsRegression(tms, adjust = 16:18)
```

Matrix B contains the 9 TMS regression coefficients, as shown in the example below:

```{r, eval = FALSE}
           bs0           bs         bi0            bi           b0
1   0.44605611 -0.172286425 -0.28845631 -0.0103719530 -0.079070482
2  -0.63479541  0.309619620  1.32187388 -0.0846722830  0.247493127
3  -0.24740092 -0.049733826  0.64390434 -0.0071768698 -0.070655046
4  -0.35290457  0.028751748  0.24985075 -0.0202048839 -0.170287969
5  -0.01846587  0.006138788 -0.23161357  0.0043496782  0.136867922
6   0.06317109  0.103245351  0.43792815 -0.0026170343  0.001647889
7   0.07581085 -0.057192403 -0.05528611 -0.0066774232  0.215741779
8  -0.01910764 -0.074350541  0.47464356 -0.0002787369 -0.154575367
9   0.13769556  0.068909513  0.83590337  0.0290674695  0.086760707
10  0.23425097 -0.050166576 -0.29430881 -0.0002431529 -0.167022626
              b1            b2         a0           a1
1  -0.0070619753 -0.0034500884 -0.9269894  0.011129428
2   0.0285864074 -0.0065340234  0.1173041  0.061815564
3   0.0082017802 -0.0001614682 -0.2955762 -0.004567065
4   0.0335833733 -0.0050643560 -1.3435844  0.008719332
5   0.0477466022 -0.0086142874 -0.4132057 -0.002285809
6   0.0201640174 -0.0005187694  2.3338259  0.089213463
7   0.0310105467 -0.0035992202 -9.8801634 -0.164766244
8  -0.0424625002  0.0058572147  0.1134623 -0.007114505
9   0.0140306025 -0.0083004128 -0.7459555  0.000597424
10  0.0008075851  0.0019799213  0.5456834 -0.007858925
```

The default TMS regression coefficients matrix in **tmsClassifier** (object `tms.coef`) is generated from 694 subjects: 273 Alzheimer's Disease patients (AD), 207 Frontotemporal Dementia patients (FTD), 67 Dementia with Lewy Bodies patients (DLB), and 147 healthy controls (HC).
