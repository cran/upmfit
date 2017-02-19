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



#' Bayesian Markov Chain Monte Carlo model building for household-community tuberculosis transmission models
#'
#' @param design.matrix A data frame type design matrix of observations for rows, variables on each observation for columns. If there are natural clusters (households), then a variable titled 'cluster' must be included. The design matrix must have exactly one column titled 'y' for outcome on each observation and at least one column of predictors.
#' @param categorical.columns A scalar (vector) specification of any column(s) of type 'factor' to be coded as a dummy variable regression. Default is 0: no columns are factor variables.
#' @param prior.beta.means A scalar (vector) specification of Gaussian prior distirbution mean(s) for the coefficient(s) in the household component of the logistic regression model. Default is noninformative on the inverse logit transformed (odds ratio) scale.
#' @param prior.beta.sd A scalar (vector) specification of Gaussian prior distirbution standard deviation(s) for the coefficient(s) in the household component of the logistic regression model.  Default is noninformative on the inverse logit transformed (odds ratio) scale.
#' @param prior.alpha A vector of prior probability parameters for the alpha term in the community component of model, a logistic regression with intercept only. Default is noninformative, with p^C being essentially uniform on the unit interval.
#' @param n.chains From JAGS documentation: number of Markov chains (default: 3).
#' @param n.iter From JAGS documentation: number of total iterations per chain (including burn in; JAGS default is 2000, updated here to 50000).
#' @param n.burnin From JAGS documentation: length of burn in, i.e. number of iterations to discard at the beginning. Default is n.iter/2, that is, discarding the first half of the simulations. If n.burnin is 0, jags() will run 100 iterations for adaption.
#' @param n.thin From JAGS documentation: thinning rate. Must be a positive integer. Set n.thin > 1 to save memory and computation time if n.iter is large. Default is max(1, floor(n.chains * (n.iter-n.burnin) / 1000)) which will only thin if there are at least 2000 simulations.
#' @return Returns a JAGS object of posterior samples for: household predictors (betas); if clusters are present, the posterior distribution of the hierarchical effects mean and standard deviation; the probability of community-acquired infection (post.comm.risk); the probability of household-acquired infection (pHHinfection); and model deviance.
#' @details The package contains long-form documentation in the form of a vignette that cover the use of the main fucntions. Use browseVignettes(package="upmfit") to access them.
#'
#'@examples
#'upmrun(design.matrix=upmdata, n.iter=100, prior.alpha = c(-2.2, 1/sqrt(3)))
#'
#' @export

upmrun<-function(design.matrix,categorical.columns=0,
                 prior.beta.means=rep(x=0,times=ncol(design.matrix[,which(names(design.matrix)!="cluster")])),
                 prior.beta.sd=rep(x=3.162278,times=ncol(design.matrix[,which(names(design.matrix)!="cluster")])),
                 prior.alpha=c(0,0.44), n.chains=3, n.iter=50000, n.burnin=n.iter/2,
                 n.thin=max(1, floor((n.iter - n.burnin) / 1000))){
   n <- nrow(suppressWarnings(upmbuilder(design.matrix, categorical.columns = categorical.columns))[[2]])
   k <- sum(suppressWarnings(upmbuilder(design.matrix, categorical.columns = categorical.columns))[[2]]["y"])
   Logit.data <- as.list(strsplit(c(names(suppressWarnings(upmbuilder(design.matrix, categorical.columns=categorical.columns))[[2]]),"n","k"), ","))
   Logit.params <- c(suppressWarnings(upmbuilder(design.matrix, categorical.columns = categorical.columns))[[3]], "post.comm.risk","pHHinfection")
   if(sum(names(suppressWarnings(upmbuilder(design.matrix, categorical.columns = categorical.columns))[[2]])=="cluster")>0) Logit.params <- c(suppressWarnings(upmbuilder(design.matrix, categorical.columns = categorical.columns))[[3]], "post.comm.risk", "pHHinfection","sigma.b.0", "mu.beta.0")
   design.matrix.rebuild<-suppressWarnings(upmbuilder(design.matrix, categorical.columns = categorical.columns))[[2]]
    for(i in 1:(length(Logit.data)-2)){
      assign(Logit.data[[i]],design.matrix.rebuild[,i])
    }
   JAGSoutput<- R2jags::jags(data=Logit.data, parameters.to.save=Logit.params, n.chains=n.chains, n.iter=n.iter, n.burnin=n.burnin, n.thin=n.thin, model.file=eval(parse(text=suppressWarnings(upmbuilder(design.matrix,categorical.columns = categorical.columns))[[1]])) )
   return(JAGSoutput)
 }

