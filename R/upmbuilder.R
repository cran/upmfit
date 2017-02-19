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
#' @param design.matrix A data frame type design matrix of observations for rows, variables on each observation for columns. If there are natural clusters (households), then a variable titled 'cluster' must be included. The design matrix must have exactly one column titled 'y' for outcome on each observation.
#' @param categorical.columns A scalar (vector) specification of any column(s) of type "factor" to be coded as a dummy variable regression. Default is 0: no columns are factor variables.
#' @param prior.beta.means A scalar (vector) specification of Gaussian prior distirbution mean(s) for the coefficient(s) in the household component of the logistic regression model. Default is noninformative on the inverse logit transformed (odds ratio) scale.
#' @param prior.beta.sd A scalar (vector) specification of Gaussian prior distirbution standard deviation(s) for the coefficient(s) in the household component of the logistic regression model.  Default is noninformative on the inverse logit transformed (odds ratio) scale.
#' @param prior.alpha A vector of prior probability parameters for the alpha term in the community component of model, a logistic regression with intercept only. Default is noninformative, with p^C being essentially uniform on the unit interval.
#' @return A list of: a JAGS model script for specification of the unified probability model to be used; the (possibly recoded) design matrix; a vector of the beta coefficient names to be used in the household component of the UPM.
#' @details The package contains long-form documentation in the form of a vignette that cover the use of the main fucntions. Use browseVignettes(package="upmfit") to access them.
#'
#'@examples
#' upmbuilder(design.matrix=upmdata)
#'
#'@examples
#' upmbuilder(design.matrix=upmdata, prior.alpha=c(-2.2,round(1/sqrt(3),3)))
#
#'@examples
#' suppressWarnings(cat(upmbuilder(design.matrix=upmdata,
#' prior.alpha=c(-2.2,round(1/sqrt(3),3)))[[1]]))
#'
#' @export


upmbuilder<-function(design.matrix, categorical.columns=0, prior.beta.means=rep(x=0,times=ncol(design.matrix[,which(names(design.matrix)!="cluster")])), prior.beta.sd=rep(x=3.162278,times=ncol(design.matrix[,which(names(design.matrix)!="cluster")])), prior.alpha=c(0,0.44) ){
    # setup -------------------------------------------------------------------
    #SET UP PARAMETERS FOR CREATION OF ADDITIONAL MATRICES, REMOVE OUTCOME Y
    if(nrow(stats::na.omit(design.matrix))!=nrow(design.matrix)){stop("Design matrix has missing values. Omit rows with missing data (using function na.omit()) before entering design matrix into function.")}
    original.matrix=design.matrix
    original.categorical.columns<-categorical.columns
    design.matrix=design.matrix[,-which(names(design.matrix)=="y")]
    if(names(original.matrix)[ncol(original.matrix)]!='y')categorical.columns=categorical.columns-1 #accounts for removal of 'y' from design matrix in error catches
    #END SET UP PARAMETERS FOR CREATION OF ADDITIONAL MATRICES

    # Error catches -------------------------------------------------------------------
    #CREATE ERROR MESSAGES:
    #if not enough columns
    if(ncol(design.matrix)<2){stop("There needs to be at least one outcome variable and one input predictor to proceed.")}

    #needs to be an outcome variable called 'y'
    if(sum(names(original.matrix)=="y")>1 | sum(names(original.matrix)=="y")==0){stop("There needs to be strictly one outcome variable titled 'y' in design matrix.")}

    #needs to be sufficient observations (rows)
    if(nrow(design.matrix)<3){stop("There needs to sufficient observations (rows) to proceed.")}

    #if there's a factor not mentioned
    for(i in 1:ncol(design.matrix)){
      if(is.factor(design.matrix[,i]) & sum(original.categorical.columns)==0 ){stop("At least one of the predictor columns in the design matrix is of type 'Factor': specify 'categorical.columns' argument.")}}

    #if variable marked as categorical in function execution is not in fact a factor, reformat:
    if(sum(categorical.columns!=0)>0 ){
      if(sum(sapply(data.frame(design.matrix[,categorical.columns]),is.factor))!=length(categorical.columns) ){
        stop("A variable(s) is listed by user as factors but data type in design matrix is not actually a factor. Reformat.")
      }
    }
    #warning on priors for categorical predictors
    if (sum(categorical.columns)>0) {
      warning("If entering non-default prior distributions on beta vector, make sure to structure the appropriate number of prior means and SDs to include multi-level categorical vectors (all non-factor covariates, plus [the total number of factor levels for each factor variable] minus [the total number of factor variables] plus the intercept), otherwise priors will be reset to noninformative.)")
    }
    #str() type error catch
    if(sum(names(original.matrix)=="cluster")>0){
      if (sum(is.numeric(original.matrix$cluster))==0) {
        stop("'cluster' variable in design matrix is not numeric. ")
      }
    }
    #warning for multiple categorical variables
    if (sum(categorical.columns==0)>1) {
      warning("Multiple categorical variables entered: intercept will be alphabetical reference group for both factors.")
    }
    #END ERROR CATCHES:

    # make categorical design matrix fomulation ----------------------------------------------
    #makes a dataset titled "design" which reformats the original "design.matrix" into the desired format for categorical variables
    #object "original.matrix" still exists, and has the original formulation of the design matrix in working memory
    if(sum(categorical.columns)>0){cats<-names(design.matrix)[categorical.columns]}

    #reformat design matrix if categorical factors
    design<-c()
    if(sum(categorical.columns)>0){
      for(i in 1:length(categorical.columns)){
        design<-cbind(design,stats::model.matrix(rep(0,nrow(design.matrix))~eval(parse(text=paste0(paste0("design.matrix$",cats[i]),collapse="+"))))[,-1])
      }
    }
    if(sum(categorical.columns)>0){
      cols<-c()
      for(j in 1:length(categorical.columns)){
        cols <- c(cols,paste(names(design.matrix)[categorical.columns][j], names(table(design.matrix[,names(design.matrix)[categorical.columns][j]]))[-1],sep=""))
      }
      colnames(design)<-cols
    }
    if(sum(categorical.columns)>0){design.matrix<-data.frame(design.matrix[,-categorical.columns],design)}
    # make categorical design matrix ------------------------------------------

    # clustering --------------------------------------------------------------
    #Define whether there is a cluster or not, and program a section in the BUGS prior probability section to deal with it if so,
    #then removes the cluster variable from the design matrix
    if(sum(names(original.matrix)=="cluster")>0){cluster<-design.matrix$cluster}
    cluster.section<-paste0("",collapse='\n')
    if(sum(names(original.matrix)=="cluster")>0) {cluster.section<-c('for(j in 1:max(cluster)){\nb.0[j]~dnorm(b.0.hat[j],tau.beta.0)\nb.0.hat[j]<-mu.beta.0\n}\nmu.beta.0~dnorm(0,0.5)\ntau.beta.0<-pow(sigma.b.0,-2)\nsigma.b.0~dunif(0,10)\n')}
    if(sum(names(original.matrix)=="cluster")>0){design.matrix<-design.matrix[,-which(names(design.matrix)=="cluster")]}
    # end clustering --------------------------------------------------------------

    # model cohesion ----------------------------------------------------------
    #define betas based on design matrix (including intercept)
    betas<-paste0('beta.', 0:(ncol(design.matrix)))

    #if priors are not user-defined, go with default, otherwise, go with defined priors THIS HAS TO BE REDONE!!
    if(sum(categorical.columns)>0 & length(prior.beta.means)<length(betas)){prior.beta.means=c(prior.beta.means[-c(which(sapply(original.matrix, is.factor))+1)], rep(x=0,times=ncol(design.matrix[,-which(names(design.matrix)%in%names(original.matrix))])-1))}
    if(sum(categorical.columns)>0 & length(prior.beta.sd)<length(betas)){prior.beta.sd=c(prior.beta.sd[-c(which(sapply(original.matrix, is.factor))+1)], rep(x=3.162278,times=ncol(design.matrix[,-which(names(design.matrix)%in%names(original.matrix))])-1))}
    b.priors<-paste0(betas,paste('~dnorm(',round(prior.beta.means,4),",",round(1/prior.beta.sd^2,4) ,')',sep=""),collapse='\n')

    #Construct model fomulation
    linear.terms <- paste0(betas[1],"+",paste0(paste0(betas[-1],'*'),colnames(design.matrix),'[i]', collapse='+'),sep="")
    if(sum(names(original.matrix)=="cluster")>0){linear.terms<-paste0(linear.terms,"+b.0[cluster[i]]")}
    formula <- noquote(paste0('logit(theta[i])<-',  linear.terms  )  )

    model.formula<-paste('function(){', 'for(i in 1:n){', 'y[i]~dbin(theta[i]+Pc[i]-Pc[i]*theta[i], 1)', formula, 'logit(Pc[i])<-alpha', '}',cluster.section, 'post.comm.risk<-exp(alpha)/(1+exp(alpha))','pHHinfection<-1-(post.comm.risk+(1-pTBinfection))' ,'pTBinfection~dbeta(1.2+k,1.2+n-k)',b.priors, paste('alpha~dnorm(',prior.alpha[1],',',prior.alpha[2],')',sep=""), '}', sep='\n')
    return(list(model.formula=model.formula, design.matrix=data.frame(design.matrix,cluster=original.matrix$cluster,y=original.matrix$y), betas=betas))
    #end model cohesion
    #end function:
  }








