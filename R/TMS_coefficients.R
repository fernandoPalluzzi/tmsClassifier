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


#' @title Default TMS regression coefficients
#'
#' @description Data.frame of TMS regression coefficients derived from 
#' default TMS data (Benussi, 2017). They include regression coefficients 
#' over time (t) for short-interval intracortical inhibition 
#' (SICI = bs0 + bs*t), intracortical facilitation (ICF = bi0 + bi*t), 
#' short-latency afferent inhibition (SAI = b0 + b1*t + b2*t^2), and 
#' long-interval intracortical inhibition (LICI = a0 + a1*t).
#' @name tms.coef
#' @usage tms.coef
#' @docType data
#' @format
#' "tms.coef" is a data.frame of 649 subjects (rows) and 9 TMS regression 
#' coefficients (columns). Subjects are subdivided into 273 Alzheimer's 
#' Disease patients (AD), 207 Frontotemporal Dementia patients (FTD), 
#' 67 Dementia with Lewy Bodies patients (DLB), and 147 healthy controls (HC).
#'
#' @references
#'
#' Benussi A, Di Lorenzo F, Dell'Era V, Cosseddu M, Alberici A, Caratozzolo S, 
#' Cotelli MS, Micheli A, Rozzini L, Depari A, Flammini A, Ponzo V, 
#' Martorana A, Caltagirone C, Padovani A, Koch G, Borroni B. (2017). 
#' Transcranial magnetic stimulation distinguishes Alzheimer disease from 
#' frontotemporal dementia. Neurology, 89(7):665-672.
#' https://doi.org/10.1212/WNL.0000000000004232
#'
NULL
