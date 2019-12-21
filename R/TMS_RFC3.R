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


#' @title AD-DLB Random Forest Classifier (RFC3)
#'
#' @description This is the third of three RFCs used by tmsClssifier to
#' predict subjects' disgnosis on the base of TMS data. For each RFC2
#' classification that yielded "Non-FTD", RFC3 decides if the respective
#' subject is of "AD" (Alzheimer's Disease) or "DLB" (Dementia with Lewy
#' Bodies) class.
#' @name TMS_RFC3
#' @usage tms.rfc3
#' @docType data
#' @format
#' "tms.rfc3" is an object of class \code{\link[randomForest]{randomForest}}.
#'
NULL
