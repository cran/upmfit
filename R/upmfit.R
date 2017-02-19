# Unified Probability Model Fitting Copyright (C) 2016-2017, Avery I. McIntosh, Boston University.
# This program is free software: you can redistribute it and/or
# modify it under the terms of the GNU General Public License as published by the
# Free Software Foundation, either version 3 of the License, or (at your option)
# any later version.  This program is distributed in the hope that it will be
# useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
# Public License for more details.  You should have received a copy of the GNU
# General Public License along with this program.  If not, see
# <http://www.gnu.org/licenses/>.

#' upmfit: Bayesian MCMC sampling with the Unified Probability Model
#'
#' The upmfit package uses methods described by McIntosh, et al. "Extensions to Bayesian Generalized Linear Mixed Effects Models for Household Tuberculosis Transmission," \emph{Statistics in Medicine} (2017), to implement functions that generate model script in the JAGS language for the Unified Probability Model (UPM), and to run a JAGS MCMC sampler for posterior probability desnity estimation with the model script.
#'
#' It should be noted that, while the UPM and R package upmfit were designed specifically with tuberculosis household contact studies in mind, the method is applicable to any scenario where there are two competing risks for a single outcome with unobserved linkage between the sources of risk and the outcome. Example situation where this method could be applicable are in profiling the risk of MRSA infection from a nosocomial versus other source, or in modeling risk of TB infection among populations utilizing homeless shelters. Non-clinical examples could be in an industrial process control setting where a manufactured item has two sources of degradation, one observed and one unobserved, or in a chemical mixing procedure where catalysis can occur from the mixing of an exogenous agent or autocatalysis.
#'
#' @section upmfit functions:
#' upmfit has: function \emph{upmbuilder()} to create the JAGS model script and display a (possibly redefined) design matrix from a design matrix input by the user; \emph{upmrun()} to generate posterior parameter estimates via a Markov Chain Monte Carlo sampler; and synthetic dataset \emph{upmdata}.
#'
#' For a longer introduction, see the introductory vignette for this package. Use command vignette("upm-primer", package="upmfit") or browseVignettes(package = "upmfit") to access the vignettes.
#'
#' @docType package
#' @name upmfit
NULL
