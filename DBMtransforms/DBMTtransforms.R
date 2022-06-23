## **************************************************************************
##
##    (c) 2018-2022 Guillaume Guénard
##        Department de sciences biologiques,
##        Université de Montréal
##        Montreal, QC, Canada
##
##    **Distance-based Mahalanobis transformations**
##
##    This file is part of codep
##
##    codep is free software: you can redistribute it and/or modify
##    it under the terms of the GNU General Public License as published by
##    the Free Software Foundation, either version 3 of the License, or
##    (at your option) any later version.
##
##    codep is distributed in the hope that it will be useful,
##    but WITHOUT ANY WARRANTY; without even the implied warranty of
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##    GNU General Public License for more details.
##
##    You should have received a copy of the GNU General Public License
##    along with codep.  If not, see <https://www.gnu.org/licenses/>.
##
##    R source code file
##
## **************************************************************************
##
#' 
#' Distance-based Mahalanobis Transformations
#' 
#' Function to obtain DNA distance-based forward and inverse linear data
#' transformation matrices.
#' 
#' @name DBMTransforms
#' 
#' @param x A DNA distance matrix (class \code{\link{dist}} object obtained
#' from package ape function \code{\link{dist.dna}}.
#' @param trf A transformation function converting the pairwise distances
#' between the DNA sequences into an estimator of the molecular covariance
#' (default: \code{function(x) -0.5*x^2}, i.e., the Gower transformation used
#' for Principal Coordinate Analysis).
#' @param tol Tolerance value for principal components rejection on the basis of
#' their absolute eigenvalue.
#' @param object A 'dna_pc' object from function \code{dna_PC}
#' @param d A vector or matrix of pairwise distances between the ESV found
#' \code{object} and a set of arbitrary objects (generally not found in
#' \code{object}).
#' @param idx Indices of the ESV used to build the DNA kernel map.
#' @param k_id Indices of the kernels chosen as representative ESV (typically,
#' medoids obtained from the DNA distance matrix).
#' @param ... Arguments to be given to the function passed as argument
#' \code{trf}.
#' 
#' @return Function \code{dna_PC} returns a \code{dna_pc}-class object
#' containing the following components:
#' \describe{
#' \item{trf}{The distance transformation function.}
#' \item{param}{A list of optional arguments for \code{trf}.}
#' \item{mn}{A list of the column means and the global mean of covariance
#' matrix}
#' \item{u}{The eigenvectors of the covariance matrix.}
#' \item{lambda}{The eigenvalues of the covariance matrix.}
#' }
#' In addition to these components, a a \code{dna_cov}-class object contains:
#' \describe{
#' \item{tr_mat}{The forward transformation matrix (from ESV to principal
#' coordinates).}
#' \item{inv_tr_mat}{The inverse transformation matrix (from principal
#' coordinates to ESV).}
#' }
#' Function \code{dna_PCscore} returns a matrix of the scores of each new ESV
#' (rows) on the \code{dna_pc}-class object's principal coordinates.
#' 
#' @details Function \code{dna_PC} calculates the principal coordinates of a
#' set of ESV from their pairwise evolutionary distances obtained the aligned
#' ESV's DNA sequences using function \code{\link{dist.dna}} from package 'ape'.
#' 
#' Given an existing \code{dna_pc}-class object built from a set of existing,
#' ESV, function \code{dna_PCscore} calculates the scores of one or more ESV
#' on the principal coordinates from the pairwise distance between the set of
#' ESV used to calculate the \code{dna_pc}-class object and the new ESV.
#' 
#' When all the non-zero eigenvalues of the principal coordinates (PC) are
#' positive, one can use function \code{dna_PCcov} to calculate the forward and
#' inverse transformation matrices associated with the implicit among-sequence
#' molecular covariance. Applied to responses associated with the ESV (such as
#' the matrix of ESV read numbers) the forward transformation matrix allows one
#' to perform the analyses (multivariate regression, redundancy analysis) on the
#' basis of the Mahalanobis metric using the implicit PC covariance as a
#' Restricted Maximum Likelihood (REML) estimate of the molecular covariance
#' among the DNA sequences. The inverse transformation matrix can be used to
#' obtain the original response values from the transformed values, possibly
#' with some information loss (depending on the number of non-zero eigenvalues,
#' which corresponds to the rank of the implicit covariance matrix).
#' 
#' @author Guillaume Guénard <guillaume.guenard@@gmail.com>
#' 
## @importFrom 
#' 
## @examples
#' 
#' @export






#' 
#' @describeIn DBMTransforms
#' 
#' DNA Principal Coordinates Scores
#' 
#' Calculate principal coordinates scores for new ESV from their DNA distances.
#' 
#' @export


#' 
#' @describeIn DBMTransforms
#' 
#' DNA Molecular Covariance Transformations
#' 
#' Calculates the forward and inverse data transformation matrices associated
#' with ESV molecular covariance.
#' 
#' @export


#' 
#' @describeIn DBMTransforms
#' 
#' DNA Distances to Kernel Coordinates
#' 
#' Calculates kernel coordinates for a set of ESV from their pairwise DNA
#' distances.
#' 
#' @export



#' 
#' @describeIn DBMTransforms
#' 
#' DNA Kernel Coordinates Scores
#' 
#' Calculate kernel coordinates scores for new ESV from their DNA distances.
#' 
#' @export




#' 
