#  tmsClassifier library
#  Copyright (C) 2019 Fernando Palluzzi; Mario Grassi
#  e-mail: <fernando.palluzzi@gmail.com>
#  University of Pavia, Department of Brain and Behavioral Sciences
#  Via Bassi 21, Pavia, 27100 Italy

#  tmsClassifier is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.

#  tmsClassifier is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.

#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.

# -------------------------------------------------------------------- #


globalVariables(c("tms.rfc1", "tms.rfc2", "tms.rfc3"))

#' S4 class to represent classifier performances.
#'
#' @slot ctable Confusion 2x2 table
#' @slot performance List of classifier performance indices
#' 
#' @importFrom methods setClass
#' 
setClass("performance", representation(ctable = "table", performance = "data.frame"))

#' S4 class to represent k-fold cross-validation results.
#'
#' @slot DvsHC Performances for Disease-Healthy classidication
#' @slot FTDvsNonFTD Performances for FTD-NonFTD classidication
#' @slot ADvsDLB Performances for AD-DLB classidication
#' 
#' @importFrom methods setClass
#' 
setClass("cross.validation", representation(DvsHC = "list", FTDvsNonFTD = "list", ADvsDLB = "list"))

#' @title Missing data imputation
#'
#' @description Impute missing Transcranial Magnetic Stimulation (TMS)
#' values, using the K-nearest neighbours (KNN) algorithm from
#' \code{\link[impute]{impute.knn}} (Troyanskaya, 2001).
#' @param tms A matrix or data.frame containing subjects as rows and TMS
#' values as columns.
#' @param sici.icf Numeric vector determining the position of
#' temporally-ordered SICI-ICF columns (SICI: short-interval intracortical
#' inhibition; ICF: intracortical facilitation). By default, they should
#' be the first 7 measures (sici.icf = 1:7; 4 for SICI and 3 for ICF),
#' taken at times (interstimulus intervals): 1, 2, 3, 5, 7, 10, 15 ms.
#' Set sici.icf to NULL to exclude these values from classification.
#' @param sai Numeric vector determining the position of
#' temporally-ordered SAI (short-latency afferent inhibition) columns.
#' By default, they should be the 4 columns following sici.icf (sai = 8:11),
#' taken at time steps (interstimulus intervals): -4, 0, 4, 8 ms. Set sai
#' to NULL to exclude these values from classification.
#' @param lici Numeric vector determining the position of temporally-ordered
#' LICI (long-interval intracortical inhibition) columns. By default, they
#' should be the 3 columns following sai (lici = 12:14), taken at time
#' steps (interstimulus intervals): 50, 100, 150 ms. Set lici to NULL to
#' exclude these values from classification.
#' @param max.na Maximum number of missing values per row (default max.na = 7).
#' @param max.r Maximum percentage of missing values per row (default max.r = 0.5).
#' @param max.c Maximum percentage of missing values per column (default max.c = 0.8).
#' @param k Number of nearest neighbours to perform imputation.
#'
#' @return The original data.frame with imputed missing values.
#'
#' @importFrom impute impute.knn
#' @export
#'
#' @references
#'
#' Troyanskaya O, Cantor M, Sherlock G, Brown P, Hastie T, Tibshirani R,
#' Botstein D, Altman RB (2001). Missing value estimation methods for
#' DNA microarrays. Bioinformatics, 17(6):520-525.
#' https://doi.org/10.1093/bioinformatics/17.6.520
#'
#' @seealso \code{\link[impute]{impute.knn}} for further details on KNN
#' imputation.
#'
tmsImpute <- function(tms, sici.icf = 1:7, sai = 8:11, lici = 12:14, max.na = 7, max.r = 0.5, max.c = 0.8, k = 10) {
	j <- c(sici.icf, sai, lici)
	Y <- tms[, j]
	n.na <- NULL
	for (i in 1:nrow(Y)) n.na[i] <- sum(is.na(Y[i,]))
	i <- which(n.na <= 7)
	Y <- as.matrix(Y[i,])
	Y <- impute.knn(Y, k = k, rowmax = max.r, colmax = max.c)$data
	#Y <- data.frame(id = tms[i, 1], as.data.frame(Y))
	return(Y)
}

#' @title Generate regression parameters from TMS data
#'
#' @description Generate regression parameters from TMS data, needed for
#' subject classification. For each subject, TMS indicators (SICI-ICF, SAI,
#' LICI) are modeled as a polynomial functions of time, in the form
#' y ~ poly(t). This function estimates two parameters for SICI
#' (SICI = bs0 + bs*t), two parameters for ICF (ICF = bi0 + bi*t), three
#' parameters for SAI (SAI = b0 + b1*t + b2*t^2), and two parameters
#' for LICI (LICI = a0 + a1*t).
#' @param tms A data.frame containing subjects as rows and TMS values as
#' columns. Optionally, the user may add covariate columns (e.g., sex, age,
#' center) to adjust TMS regression estimates.
#' @param sici.icf Numeric vector determining the position of
#' temporally-ordered SICI-ICF columns (SICI: short-interval intracortical
#' inhibition; ICF: intracortical facilitation). By default, they should
#' be the first 7 measures (sici.icf = 1:7; 4 for SICI and 3 for ICF),
#' taken at times (interstimulus intervals): 1, 2, 3, 5, 7, 10, 15 ms.
#' Set sici.icf to NULL to exclude these values from classification.
#' @param sai Numeric vector determining the position of
#' temporally-ordered SAI (short-latency afferent inhibition) columns.
#' By default, they should be the 4 columns following sici.icf (sai = 8:11),
#' taken at time steps (interstimulus intervals): -4, 0, 4, 8 ms. Set sai
#' to NULL to exclude these values from classification.
#' @param lici Numeric vector determining the position of temporally-ordered
#' LICI (long-interval intracortical inhibition) columns. By default, they
#' should be the 3 columns following sai (lici = 12:14), taken at time
#' steps (interstimulus intervals): 50, 100, 150 ms. Set lici to NULL to
#' exclude these values from classification.
#' @param adjust Numeric vector determining the position of covariates to
#' adjust for. By default, adjust = NULL (no covariate adjustment is done).
#'
#' @return A data.frame of estimated regression parameters.
#'
#' @importFrom stats lm
#' @export
#'
tmsRegression <- function(tms, sici.icf = 1:7, sai = 8:11, lici = 12:14, adjust = NULL) {
	if (!is.null(sici.icf)) {
		Y <- as.matrix(tms[, sici.icf])
		t <- c(1, 2, 3, 5, 7, 10, 15)
		# SICI
		cat("Parameters estimation for SICI-ICF\n")
		bs0 <- vector()
		bs <- vector()
		for (i in 1:nrow(Y)) {
			BLR1 <- lm(Y[i, 1:4] ~ t[1:4])$coefficients
			bs0[i] <- BLR1[1]
			bs[i] <- BLR1[2]
		}
		BLR <- cbind(bs0, bs)
		# ICF
		bi0 <- vector()
		bi <- vector()
		for (i in 1:nrow(Y)) {
			BLR2 <- lm(Y[i, 5:7] ~ t[5:7])$coefficients
			bi0[i] <- BLR2[1]
			bi[i] <- BLR2[2]
		}
		# SICI.ICF
		BLR <- cbind(BLR, bi0, bi)
	}
	if (!is.null(sai)) {
		# SAI
		cat("Parameters estimation for SAI\n")
		Y <- as.matrix(tms[, sai])
		t <- c(-4, 0, 4, 8)
		b0 <- vector()
		b1 <- vector()
		b2 <- vector()
		for (i in 1:nrow(Y)) {
			QR <- lm(Y[i,] ~ t + I(t^2))$coefficients
			b0[i] <- QR[1]
			b1[i] <- QR[2]
			b2[i] <- QR[3]
		}
		QR <- cbind(b0, b1, b2)
	}
	if (!is.null(lici)) {
		# LICI
		cat("Parameters estimation for LICI\n")
		Y <- as.matrix(tms[, lici])
		t <- c(50, 100, 150)
		a0 <- vector()
		a1 <- vector()
		for (i in 1:nrow(Y)) {
			LR <- lm(Y[i,] ~ t)$coefficients
			a0[i] <- LR[1]
			a1[i] <- LR[2]
		}
		LR <- cbind(a0, a1)
	}
	B <- cbind(BLR, QR, LR)

	# Covariates adjustment
	if (is.null(adjust)) {
		cat("Done.\n")
		return(B)
	} else {
		cat("Covariates adjustment\n")
		L <- length(adjust)
		Z <- tms[, adjust]
		X <- "model.matrix(~ Z[,1]"
		if (L > 1) {
			for (i in 2:L) {
				X <- paste(X, " + Z[,", i, "]", sep = "")
			}
		}
		X <- eval(parse(text = paste(X, ")", sep = "")))
		Q <- diag(nrow(X)) - X %*% solve(t(X) %*% X) %*% t(X)
		R <- Q %*% B
		cat("Done.\n")
		return(R)
	}
}

#' @title TMS Random Forest Classifier (RFC)
#'
#' @description Classify subjects with unknown diagnosis in one class among
#' Healthy subject (HC), Frontotemporal Dementia (FTD), Alzheimer's Disease
#' (AD), or Dementia with Lewy Bodies (DLB), using estimated TMS regression
#' parameters. This function provides a default RFC, requiring an input
#' data.frame of 9 columns (4 for SICI-ICF, 3 for SAI, and 2 for LICI; see
#' \code{\link[tmsClassifier]{tmsRegression}} for details). See function
#' \code{\link[tmsClassifier]{buildPredictor}} to build a custom RFC.
#' Alternatively, if a diagnosis column is present in the input data.frame,
#' it can be used to test RFC performances through cross-validation.
#' @param data A data.frame containing subjects as rows and TMS parameters
#' as columns. Optionally, a column specifying patient diagnosis can be
#' added to test RFC performances using cross-validation.
#' @param cv A logical value. If cv is TRUE, the cross-validation mode is
#' enabled, and a status (i.e., diagnosis) column is needed (see below).
#' By default, cv = FALSE and the diagnosis prediction mode is enabled.
#' @param status Numeric value indicating the column for subject diagnosis
#' (default status = 10). This value is ignored if cv is set to FALSE.
#' @param p Proportion of the input dataset to be used for the
#' cross-validation. For instance, if p = 0.75, a randomly chosen 25% of
#' the input data.frame will be isolated and returned as validation set,
#' while the remaining 75% is used for the cross-validation. Every partition
#' is balanced so that the original class proportion is preserved. By default,
#' p = 1 (i.e., the entire dataset is used for cross-validation) and no
#' validation set is generated.
#' @param k Number of cross-validation iterations (default k = 5).
#'
#' @return If cv = FALSE, the vector of diagnosis prediction for each
#' subject is returned. If cv = TRUE, an object of class cross.validation
#' is generated, including cross-validation results. If p < 1, a list of
#' three objects is returned:
#' \enumerate{
#' \item "CV", cross-validation performances;
#' \item "trainingSet", dataset partition used for the cross-vaidation;
#' \item "validationSet", dataset partition isolated before cross-validation.
#' }
#'
#' @importFrom CMA GenerateLearningsets classification
#' @importFrom methods new
#' @export
#'
#' @references
#'
#' Slawski M, Daumer M, Boulesteix AL (2008). CMA - a comprehensive 
#' Bioconductor package for supervised classification with high 
#' dimensional data. BMC Bioinformatics, 9:439. 
#' https://doi.org/10.1186/1471-2105-9-439
#'
#' @seealso \code{\link[CMA]{GenerateLearningsets}} and
#' \code{\link[CMA]{classification}}.
#'
tmsClassify <- function(data, cv = FALSE, status = 10, p = 1, k = 5) {

	if (cv) {

		# Data preparation

		cat("\nData preparation ...\n")
		data <- recode(data = data, status = status)

		if (p < 1) {
			V <- vpart(data, p = p, status = status)
			data <- V$training.set
			V <- V$validation.set
		}

		# Diseased vs. HC

		cat("\nDiseased vs Healthy Controls classification\n")
		#set.seed(123)
		LS1 <- GenerateLearningsets(n = nrow(data), data$D_HC,
									method = "CV",
									fold = k,
									strat = TRUE)
		RF1 <- classification(as.matrix(data[, 1:(status-1)]), data$D_HC,
							  learningsets = LS1,
							  classifier = "rfCMA")

		C1 <- list()
		for (i in 1:length(RF1)) {
			C1[[i]] <- performance(RF1[[i]]@y, RF1[[i]]@yhat)
		}

		# FTD vs. NonFTD

		cat("\nFTD vs non-FTD classification\n")
		X2 <- data[data[, status] != "HC",]

		#set.seed(123)
		LS2 <- GenerateLearningsets(n = nrow(X2), X2$FTD_NonFTD,
									method = "CV",
									fold = k,
									strat = TRUE)
		RF2 <- classification(as.matrix(X2[, 1:(status-1)]), X2$FTD_NonFTD,
							  learningsets = LS2,
							  classifier = "rfCMA")

		C2 <- list()
		for (i in 1:length(RF2)) {
			C2[[i]] <- performance(RF2[[i]]@y, RF2[[i]]@yhat)
		}

		# AD vs. DLB

		cat("\nAD vs DLB classification\n")
		X3 <- data[data[, status] %in% c("AD", "DLB"),]

		#set.seed(123)
		LS3 <- GenerateLearningsets(n = nrow(X3), X3$AD_DLB,
									method = "CV",
									fold = k,
									strat = TRUE)
		RF3 <- classification(as.matrix(X3[, 1:(status-1)]), X3$AD_DLB,
							  learningsets = LS3,
							  classifier = "rfCMA")

		C3 <- list()
		for (i in 1:length(RF3)) {
			C3[[i]] <- performance(RF3[[i]]@y, RF3[[i]]@yhat)
		}

		CV <- new("cross.validation", DvsHC = C1, FTDvsNonFTD = C2, ADvsDLB = C3)
		cat("\nDone.\n")
		if (p < 1) {
			return(list(crossValidation = CV, trainingSet = data, validationSet = V))
		} else {
			return(CV)
		}

	} else {

		n <- nrow(data)
		v <- vector()
		for (i in 1:n) {
			cat("\nSubject", i, "of", n, "\n  ")
			v <- c(v, predict.subject(data[i,]))
			cat("\n")
		}
		return(v)
	}
}

recode <- function(data, status) {

	data$D_HC <- "Diseased"
	data$D_HC[data[, status] == "HC"] <- "HC"
	data$D_HC <- as.factor(data$D_HC)

	data$FTD_NonFTD <- "NonFTD"
	data$FTD_NonFTD[data[, status] == "FTD"] <- "FTD"
	data$FTD_NonFTD <- as.factor(data$FTD_NonFTD)

	data$AD_DLB <- "AD"
	data$AD_DLB[data[, status] == "DLB"] <- "DLB"
	data$AD_DLB <- as.factor(data$AD_DLB)

	return(data)
}

#' @title Remove outliers using subject-wise Brier scores
#'
#' @description Detect and remove outliers from a dataset using per subject
#' Brier scores (Brier, 1950).
#' @param data A matrix or data.frame containing subjects as rows and TMS
#' parameters as columns. A further column specifying subject diagnosis
#' is required to run RFC cross-validation.
#' @param status Numeric value specifying subject diagnosis (default
#' status = 10).
#' @param k Number of cross-validation cycles (default k = 5).
#' @param b Maximum Brier score beyond which a subject is considered as
#' an outlier. By default, b = 1. If b is set to NULL, outliers will be
#' automatically set to Q3 + 1.5*(Q3 - Q1).
#'
#' @return The input data.frame without outliers.
#'
#' @importFrom CMA GenerateLearningsets classification evaluation
#' @importFrom graphics boxplot
#' @export
#'
#' @references
#'
#' Brier GW (1950). Verification of forecasts expressed in terms of probability.
#' Monthly Weather Review, 78(1):1-3.
#' https://doi.org/10.1175/1520-0493(1950)078<0001:VOFEIT>2.0.CO;2
#'
#' @seealso \code{\link[CMA]{evaluation}}
#'
rmout <- function(data, status = 10, k = 5, b = 1) {

	L <- GenerateLearningsets(n = nrow(data), data[, status],
	                          method = "CV",
	                          fold = k,
	                          strat = TRUE)
	C <- classification(as.matrix(data[, 1:(status-1)]), data[, status],
	                    learningsets = L,
	                    classifier = "rfCMA")

	E <- evaluation(C, scheme = "observationwise", measure = "brier score")

	if (is.null(b)) {
		outliers <- boxplot(E@score, plot = FALSE)$out
		return(data[-which(E@score >= min(outliers)),])
	} else {
		return(data[-which(E@score > b),])
	}
}

#' @title Build a RFC and test its performances
#'
#' @description Build a RFC from TMS regression parameters and test its
#' classification performances.
#' @param data A matrix or data.frame containing subjects as rows and TMS
#' parameters as columns. A further column specifying subject diagnosis
#' can be given to validate the classifier.
#' @param status Numeric value specifying subject diagnosis (default
#' status = 10).
#' @param vset An optional data.frame of the same format of the input data.
#' This will be used as external validation set. The validation set can
#' be generated from input data, using \code{\link[tmsClassifier]{tmsClassify}}.
#' @param n Number of bootstrap sampled trees of the RFC (default n = 10000).
#' Increasing n improves RFC classification performance, but will sensibly
#' increase computational time.
#' @param m Number of random variables for each tree split (default m = 3).
#' The suggested value is the squared root of the total number of variables
#' (the default value is set to 3 since 9 TMS parameter are expected by
#' the default classifier).
#'
#' @return An object of class \code{\link[randomForest]{randomForest}}.
#' An optional object of class \code{\link[tmsClassifier]{performance}}
#' is returned, if vset is not NULL.
#'
#' @importFrom randomForest randomForest
#' @importFrom stats predict
#' @export
#'
#' @references
#'
#' Liaw A, Wiener M. Classification and Regression by randomForest (2002). 
#' R News, 2(3):18-22. https://doi.org/10.1023/A:1010933404324
#'
#' @seealso \code{\link[randomForest]{randomForest}}
#'
buildPredictor <- function(data, status = 10, vset = NULL, n = 10000, m = 3) {
	x <- colnames(data)
	model <- paste("as.factor(", x[status], ") ~ ", paste(x[-status:-length(x)],
	               collapse = " + "), collapse = "", sep = "")

	RFC <- randomForest(eval(parse(text = model)), data = data,
						ntree = n,
						mtry = m,
						keep.forest = TRUE,
						proximity = TRUE,
						importance = TRUE)

	if (is.null(vset)) {
		return(RFC)
	} else {
		C <- predict(RFC, vset)
		P <- performance(vset[, status], C)
		return(list(RFC = RFC, performance = P))
	}
}

#' @title Data partitioning utility
#'
#' @description Extract a random partition from an input dataset.
#' @param data An input data.frame.
#' @param p Proportion of input rows to be extracted from the input
#' dataset.
#' @param status Numeric value indicating the column for subject diagnosis
#' (default status = NULL). Allowed values for this column are "HC", "FTD",
#' "AD", and "DLB". If status is not NULL, random partitioning will preserve
#' the status column proportions.
#' @param shuffle A logical value. If TRUE, the input rows are randomly
#' shuffled before data partitioning.
#'
#' @return A list of 2 data.frames:
#' \enumerate{
#' \item "training.set", the portion of the input data defined by p;
#' \item "validation.set", the portion of the input data defined by 1-p.
#' }
#'
#' @export
#'
vpart <- function(data, p, status = NULL, shuffle = FALSE) {
	if (shuffle) {
		data <- data[sample(1:nrow(data), replace = FALSE),]
	}
	if (!is.null(status)) {
		HC <- data[data[, status] == "HC",]
		FTD <- data[data[, status] == "FTD",]
		AD <- data[data[, status] == "AD",]
		DLB <- data[data[, status] == "DLB",]
		h <- sample(nrow(HC), size = p*nrow(HC), replace = FALSE)
		f <- sample(nrow(FTD), size = p*nrow(FTD), replace = FALSE)
		a <- sample(nrow(AD), size = p*nrow(AD), replace = FALSE)
		d <- sample(nrow(DLB), size = p*nrow(DLB), replace = FALSE)
		X <- rbind(HC[h,], FTD[f,], AD[a,], DLB[d,])
		V <- rbind(HC[-h,], FTD[-f,], AD[-a,], DLB[-d,])
	} else {
		j <- sample(nrow(data), size = p*nrow(data), replace = FALSE)
		X <- data[j,]
		V <- data[-j,]
	}
	return(list(training.set = X, validation.set = V))
}

predict.subject <- function(data) {
	C1 <- predict(tms.rfc1, data)
	if (C1 == 0) {
		cat("Healthy subject -> HC\n")
		return("HC")
	} else {
		x <- "Diseased subject"
		C2 <- predict(tms.rfc2, data)
		if (C2 == 1) {
			cat(paste(x, "FTD\n", sep = " -> "))
			return("FTD")
		} else {
			C3 <- predict(tms.rfc3, data)
			if (C3 == 1) {
				cat(paste(x, "AD\n", sep = " -> "))
				return("AD")
			} else {
				cat(paste(x, "DLB\n", sep = " -> "))
				return("DLB")
			}
		}
	}
}

#' @title Evaluate classifier performances
#'
#' @description Compute performance indices for a binary classifier, from
#' a vector of observed values and a vector of predictions.
#' @param obs Vector of observed values.
#' @param pred Vector of predicted values.
#' @param CT An object of class "table". Alternatively to the vectors of
#' observed and predicted values, the user could provide a 2x2 contingency
#' table.
#' @param y Outcome encoding. Use "0,1" ("1,2") if the binary encoding is
#' 0 (1) for controls and 1 (2) for cases.
#'
#' @return An object of class \code{\link[tmsClassifier]{performance}}.
#'
#' @importFrom methods setClass
#' @export
#'
performance <- function(obs = NULL, pred = NULL, CT = NULL, y = "0,1")
{
	# Confusion Table
	if (class(CT) != "table") {
		CT <- table(obs, pred)
	}

	#print(CT)
	#cat("\n")
	if (y == "1,2") {
		a <- CT[1, 1]
		b <- CT[1, 2]
		c <- CT[2, 1]
		d <- CT[2, 2]
	}
	if (y == "0,1") {
		a <- CT[2, 2]
		b <- CT[2, 1]
		c <- CT[1, 2]
		d <- CT[1, 1]
	}
	Se <- a/(a + b)                  # Sensitivity, TP/(TP + FN)
	Sp <- d/(c + d)                  # Specificity, TN/(TN + FP)
	Ac <- (a + d)/(a + b + c + d)    # Accuracy
	R <- a/(a + b)                   # Recall, TP/(TP + FN)
	P <- a/(a + c)                   # Precision, TP/(TP + FP)
	H <- (2*R*P)/(R + P)             # Harmonic accuracy (F1-measure)
	G <- sqrt(R*P)                   # Geometric accuracy

	P <- new("performance", ctable = CT, performance = data.frame(Se, Sp, Ac, P, H, G))

	return(P)
}

