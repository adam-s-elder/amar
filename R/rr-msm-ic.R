#' Estimate both the parameter, and the influence
#' curves used for estimating the projected risk ratio from a working
#' marginal structural mean model.
#'
#' The function takes a dataset, and will returns the estimate for the
#' parameter, and or the estimated influence curve at each observation.
#' The first column of obs_data should be the binary outcome of interest.
#'
#' Most of this code has been copied and slightly modified from the
#' ltmle package and is copyrighted by Joshua Schwab under the terms of
#' the GPL-2 license.
#'
#' Obtain an estimator of the probability delta = 1 given w
#' @param obs_data the observed data.  The first column should be the outcome.
#' @param what the desired return value. Should be one of `"ic"`
#' (influence curve), `"est"` (estimate), or `"both"`.
#' @param control any other control parameters to be passed to the estimator.
#' @return If `what` is
#'
#' - `"est"`, then return the estimated parameter.
#'
#' - `"ic"`, then return the estimated IC of the parameter estimate.
#'
#' - `"both"`, then return both the parameter estimate and
#' corresponding estimated IC.
#' @examples
#' #not run (make sure to load in SuperLearner if running)
#' #set.seed(1010)
#' #fake_dat <- data.frame(y = rbinom(100, size = 1, prob = 0.5),
#' #                      a = rbinom(100, size =  1, prob = 0.5),
#' #                       w = matrix(rnorm(500), ncol = 5))
#' # rr.msm.ic(fake_dat)
#'
#' @export

rr.msm.ic <- function(obs_data, what = "both", control = NULL){
  if (!(what %in% c("ic", "est", "both"))) {
    stop("what must be one of ic (influence curve), est (estimate), or both")
  }
  ret <- list()
  w_covs <- which(!colnames(obs_data) %in% c("y", "a"))
  fin_IC <- matrix(NA, nrow = nrow(obs_data),
                   ncol = ncol(obs_data) - 2)
  psi_hat <- rep(NA, length(w_covs))
  w_covs <- which(!colnames(obs_data) %in% c("y", "a"))
  fin_IC <- matrix(NA, nrow = nrow(obs_data),
                   ncol = ncol(obs_data) - 2)
  for (cov_idx in 1:length(w_covs)) {
    cov_IC <- my_tmleMSM(Y = obs_data$y,
                         V = obs_data[, w_covs[cov_idx]],
                         A = obs_data$a,
                         W = obs_data[, w_covs[-cov_idx]],
                         MSM = "A*V", gform = A~1,
                         g.SL.library = NULL,
                         Q.SL.library = list(#"SL.glm",
                           c("SL.glmnet", "screen.glmnet")),
                         hAVform = A~1, family = "binomial",
                         ret_IC = TRUE)
    fin_IC[, cov_idx] <- cov_IC$IC[, 4]
    psi_hat[cov_idx] <- cov_IC$psi[4]
  }
  if (what %in% c("both", "est")) {
    ret$est <- psi_hat
  }
  if (what %in% c("both", "ic")) {
    ret$ic <- fin_IC
  }
  return(ret)
}

#----- function estQcvSL ----
# purpose: Obtain cross-validated estimates for initial Q using discrete SL.
# 	This function will return cross-validated predicted initial values
#   corresponding to the best algorithm in the library, or the convex combination
#   calculated by SL itself, as chosen by cvRSS.
#   The fitted value for each observation i, is predicted from a fit based on a training set
#	excluding the fold containing observation i. Observations with same id are grouped in the
#   same fold.
# arguments:
# 	Y - outcome
# 	X - design matrix with (Z,A,W)
# 	SL.library - prediction algorithms
# 	V - number of outer cv-folds for disscrete SL
# 	V_SL - number of folds for internal SL cross-validation
# 	family - binomial or gaussian
#   Delta - missingness indicator
#   Qbounds - bounds on predicted values for Q
#   id - subject identifier
# returns:
#  Q - nx5 matrix of predicted values on linear scale (e.g. logit if Qfamily=binomial)
#-----------------------------------------------

estQcvSL <- function(Y,X,SL.library=NULL, V=5, V_SL=5, family="gaussian", Delta, Qbounds, id, verbose){
  SL.version <- 2
  Q <- cvRSS <- best_alg <- NULL
  n <- length(Y)
  u.id <- unique(id)
  n.id <- length(u.id)
  fold <- by(sample(1:n.id),rep(1:V, length.out=n.id),function(x){which(id %in% u.id[x])})
  n_predictors <-length(SL.library)
  CDE <- length(unique(X[,1])) > 1

  if(NCOL(X) > 1){
    newX <- rbind(X,X,X,X,X)
    newX[(n+1):(3*n),1]   <- 0
    newX[(3*n+1):(5*n),1] <- 1
    newX[(n+1):(2*n),2]   <- 0
    newX[(2*n+1):(3*n),2] <- 1
    newX[(3*n+1):(4*n),2] <- 0
    newX[(4*n+1):(5*n),2] <- 1
    # We'll create a matrix of predictions - one column for each predictor in the library
    # plus one more for SL itself, with 5*n predicted values per column, corresponding to newX.
    predictions <- matrix(nrow=5*n, ncol=length(SL.library)+1)
    m_SL <- NULL
    for (v in 1:V) {
      fold_rows <- c(fold[[v]], fold[[v]]+n, fold[[v]]+2*n, fold[[v]]+3*n, fold[[v]]+4*n)
      if (class(m_SL) != "try-error") {
        train.observed <- (1:n)[-fold[[v]]][Delta[-fold[[v]]]==1]
        if(utils::packageDescription("SuperLearner")$Version < SL.version){
          arglist <- list(Y=Y[train.observed], X=X[train.observed,], newX=newX[fold_rows,],
                          V=V_SL, save.fit.library=FALSE, family=family,SL.library=SL.library,id=id[train.observed])
        } else {
          arglist <- list(Y=Y[train.observed], X=X[train.observed,], newX=newX[fold_rows,],
                          cvControl=list(V=V_SL), control = list(saveFitLibrary=FALSE), family=family,SL.library=SL.library,id=id[train.observed])
        }

        suppressWarnings(
          m_SL <- try(do.call(SuperLearner::SuperLearner, arglist))
        )
      }
      if (class(m_SL) != "try-error"){
        predictions[fold_rows,1] <- m_SL$SL.predict
        for (s in 1:n_predictors){
          predictions[fold_rows,s+1] <- m_SL$library.predict[,s]
        }
        predictions <- .bound(predictions, Qbounds)
      } else {
        stop("Super Learner failed when estimating Q. Exiting program\n")
      }
    }
    cvRSS <- colSums(Delta*(Y-predictions[1:n,])^2)
    names(cvRSS) <- c("SL", SL.library)
    best <- which.min(cvRSS)
    best_alg <- c("SL", SL.library)[best]
    Q <- matrix(data=predictions[,best], nrow=n, ncol=5, byrow=FALSE)
    colnames(Q) <- c("QAW", "Q0W", "Q1W", "Q0W.Z1", "Q1W.Z1")
  }
  if(verbose){cat("\tDiscrete SL: best algorithm = ", best_alg,"\n")}
  if (is.null(Q) | inherits(m_SL, "try-error")){
    Q <- 0
    class(Q) <- "try-error"
  }
  Qinit <- list(Q=Q, family=family, SL.library=SL.library, cvRSS=cvRSS, best_alg=best_alg)
  return(Qinit)
}

## Needed Functions:


my_tmleMSM <- function(Y,A,W,V,T=rep(1,length(Y)), Delta=rep(1, length(Y)), MSM, v=NULL,
                       Q=NULL, Qform=NULL,
                       Qbounds=c(-Inf, Inf), Q.SL.library=c("SL.glm", "SL.glmnet"),
                       cvQinit = FALSE,
                       hAV=NULL, hAVform=NULL,
                       g1W = NULL, gform=NULL,
                       pDelta1=NULL, g.Deltaform=NULL,
                       g.SL.library = NULL, # c("SL.glm", "SL.step", "SL.glm.interaction")
                       ub = 1/0.025,
                       family="gaussian", fluctuation="logistic", alpha  = 0.995,
                       id=1:length(Y), V_SL=5, inference=TRUE, verbose=FALSE, ret_IC = TRUE) {
  Y[is.na(Y)] <- 0
  n <- length(Y)
  n.id <- length(unique(id))
  if(is.null(v)){
    I.V <- rep(1, length(Y))
  } else {
    I.V <- as.numeric(V==v)
  }
  V <- as.matrix(V)
  colnames(V) <- .setColnames(colnames(V), NCOL(V), "V")

  if (sum(!sapply(W, is.numeric)) > 0) {
    stop("Currently, only numeric variables are allowed.  Please convert any character or factor variables to numeric.")
  }

  W <- as.matrix(W)
  colnames(W) <- .setColnames(colnames(W), NCOL(W), "W")

  if(identical (family, stats::binomial)){
    family <- "binomial"
  } else if (identical(family, stats::gaussian)){
    family <- "gaussian"
  }
  if(!.verifyArgs(Y,Z=NULL,A,cbind(V,W,T),Delta, Qform, gform, hAVform, g.Deltaform)){
    stop()
  }
  maptoYstar <- fluctuation=="logistic" | family=="binomial"

  #---- Stage 1 -----
  stage1 <- .initStage1(Y=Y,A=A, Q=Q, Delta=Delta, Qbounds=Qbounds, alpha=alpha, maptoYstar=maptoYstar, family=family)
  Qinit <- suppressWarnings(estimateQ(Y=stage1$Ystar,Z=rep(1, length(Y)), A=A,
                                      W=cbind(W,V,T), Delta=(I.V==1 & Delta==1),
                                      Q=stage1$Q, Qbounds=stage1$Qbounds, Qform=Qform, maptoYstar = maptoYstar,
                                      SL.library=Q.SL.library, cvQinit=cvQinit, family=family, id=id, V = V_SL, verbose=verbose))

  #---- Stage 2 -----
  if(is.null(hAV)){
    gAV <- suppressWarnings(estimateG(d=data.frame(A,V,T), hAV, hAVform,g.SL.library, id, V=V_SL, verbose,
                                      message="h(A,V)", outcome="A"))
    hAV <- cbind((1-A)*(1-gAV$g1W) + A*gAV$g1W, 1-gAV$g1W, gAV$g1W)
  } else {
    hAV <- cbind((1-A)*(hAV[,1]) + A*hAV[,2], hAV)
    gAV <- NULL
    gAV$g1W <- hAV
    gAV$type <- "User-supplied values"
    gAV$coef <- NA
  }
  colnames(hAV) <- c("hAV", "h0V", "h1V")
  if (is.null(v)){
    g <- suppressWarnings(estimateG(d=data.frame(A,V,W,T), g1W, gform,g.SL.library, id, V=V_SL, verbose,
                                    message="treatment mechanism", outcome="A"))
  } else {
    g <- suppressWarnings(estimateG(d=data.frame(A,V,W,T), g1W, gform,g.SL.library, id, V=V_SL, verbose,
                                    message="treatment mechanism", outcome="A", newdata=data.frame(A,V=v, W,T)))
  }
  g$bound <- c(0,ub)
  if(g$type=="try-error"){
    stop("Error estimating treatment mechanism (hint: only numeric variables are allowed)")
  }
  g.Delta <- estimateG(d=data.frame(Delta, Z=1, A, W,V,T), pDelta1, g.Deltaform,
                       g.SL.library,id=id, V = V_SL, verbose, "missingness mechanism", outcome="D")
  g1VW <- g$g1W * g.Delta$g1W[,"Z0A1"]
  g0VW <- (1-g$g1W) * g.Delta$g1W[,"Z0A0"]
  gAVW <- A*g1VW + (1-A)*g0VW

  MSMformula <- stats::formula(paste("Y~",MSM))

  mfA <- stats::model.frame(MSMformula, data=data.frame(Y,A,V,W,T))  # Added Y Oct 24 ,2011
  mf0 <- stats::model.frame(MSMformula, data=data.frame(Y,A=rep(0,n),V,W,T))
  mf1 <- stats::model.frame(MSMformula, data=data.frame(Y,A=rep(1,n),V,W,T))

  covar.MSM   <- stats::model.matrix(MSMformula, mfA)
  covar.MSMA0 <- stats::model.matrix(MSMformula, mf0)
  covar.MSMA1 <- stats::model.matrix(MSMformula, mf1)

  if(verbose){cat("\tTargeting Q\n")}
  C1 <- I.V * .bound(hAV[,"hAV"]/gAVW, c(0,ub)) * covar.MSM
  suppressWarnings(
    epsilon <- stats::coef(stats::glm(stage1$Ystar ~ -1 + offset(Qinit$Q[,"QAW"]) + C1, subset=(I.V==1 & Delta==1), family=Qinit$family))
  )
  Qstar <- cbind(Qinit$Q[,"QAW"] + C1 %*% epsilon,
                 Qinit$Q[,"Q0W"] + .bound(hAV[,"h0V"]/g0VW, c(0,ub)) * covar.MSMA0 %*% epsilon,
                 Qinit$Q[,"Q1W"] + .bound(hAV[,"h1V"]/g1VW, c(0,ub)) * covar.MSMA1 %*% epsilon)
  if(identical(Qinit$family, "binomial")){
    Qstar <- stats::plogis(Qstar)*diff(stage1$ab)+stage1$ab[1]
    Qinit$Q <- stats::plogis(Qinit$Q)*diff(stage1$ab)+stage1$ab[1]
    Ystar <- stage1$Ystar*diff(stage1$ab)+stage1$ab[1]
  }
  colnames(Qstar) <- c("QAW", "Q0W", "Q1W")

  if(verbose){cat("\tEvaluating MSM parameters\n")}
  d.Qstar <- data.frame(Y=c(Qstar[,"Q0W"], Qstar[,"Q1W"]),
                        rbind(mf0, mf1),
                        wts=c(hAV[,"h0V"], hAV[,"h1V"]))
  suppressWarnings(
    psi.Qstar <- stats::coef(stats::glm(MSMformula, data=d.Qstar,
                                        weights=d.Qstar$wts, family=family))
  )
  d.Qinit <- replace(d.Qstar,1, c(Qinit$Q[,"Q0W"], Qinit$Q[,"Q1W"]))
  suppressWarnings(
    psi.Qinit <- stats::coef(stats::glm(MSMformula, data=d.Qinit, weights=d.Qinit$wts, family=family))
  )

  if(inference){
    if(verbose){cat("\tCalculating variance-covariance matrix\n")}
    if(family=="binomial"){
      mAV <- stats::plogis(cbind(covar.MSMA0 %*% psi.Qstar, covar.MSMA1 %*% psi.Qstar))
    } else {
      mAV <- cbind(covar.MSMA0 %*% psi.Qstar, covar.MSMA1 %*% psi.Qstar)
    }
    colnames(mAV) <- c("m0V", "m1V")

    sigma <- calcSigma(hAV, gAVW, Y, Qstar, mAV, covar.MSM, covar.MSMA0, covar.MSMA1, I.V, Delta, ub, id, family,
                       return_IC = ret_IC)
    if(ret_IC){
      IC <- sigma$IC
      sigma <- sigma$sigma/n.id
    }else{
      sigma <- sigma/n.id
      IC <- NULL
    }
    se <- sqrt(diag(sigma))
    pvalue <- 2*stats::pnorm(-abs(psi.Qstar/se))
    lb <- psi.Qstar - 1.96 * se
    ub <- psi.Qstar + 1.96 * se
  } else {
    sigma <- se <- lb <- ub <- pvalue <- IC <- NULL
  }
  Qinit$Q <- Qinit$Q[,-1]
  returnVal <- list(psi=psi.Qstar, sigma=sigma,se=se, pvalue=pvalue, lb=lb, ub=ub, epsilon=epsilon,  psi.Qinit=psi.Qinit,  Qstar=Qstar[,-1], Qinit=Qinit, g=g, g.AV=gAV, g.Delta=g.Delta, IC=IC)
  class(returnVal) <- "tmleMSM"
  return(returnVal)
}



#---------- function .setColnames ---------------
# assign names to every unnamed column of x
# arguments
# 	x.colnames - current column names
#	x.ncols - current number of columns
# 	firstChar - prefix for internally assigned name
# return the names
#-----------------------------------------
.setColnames <- function(x.colnames, x.ncols, firstChar){
  if(is.null(x.colnames)) {
    if(x.ncols > 1){
      x.colnames <- paste(firstChar,1:x.ncols, sep="")
    } else {
      x.colnames <- firstChar
    }
  } else {
    invalid.name <- nchar(x.colnames) == 0
    if(any(invalid.name)){
      x.colnames[invalid.name] <- paste(".internal",firstChar, which(invalid.name), sep="")
    }
  }
  return(x.colnames)
}

#-------------.verifyArgs------------------
# initial checks on data passed in
#-------------------------------------------
.verifyArgs <- function(Y,Z,A,W,Delta, Qform, gform, g.Zform,g.Deltaform){
  formulas <- list(Qform, gform, g.Zform, g.Deltaform)
  validFormula <- sapply(formulas, function(x){identical(class(try(stats::as.formula(x))), "formula")})
  validNames <- c("Y", "Z", "A", ".", "Delta", colnames(W))
  validTerms <- rep(TRUE, length(formulas))
  validTerms[validFormula] <- sapply(formulas[which(validFormula)],
                                     function(x){is.null(x) || all(all.names(stats::as.formula(x), functions=FALSE) %in% validNames)})
  ok <- c(length(Y) == length(A) & length(A) == NROW(W) & length(A)==NROW(Delta),
          sum(is.na(A), is.na(W)) == 0,
          all(A[!is.na(A)] %in% 0:1),
          is.null(Z) || all(Z[!is.na(Z)] %in% 0:1),
          is.null(Z) || length(unique(Z)) == 1 | (length(unique(Z)) > 1 & length(unique(A))>1),
          validFormula,
          validTerms
  )
  warning_messages <- c("\t'Y', 'A', 'W', 'Delta', must contain the same number of observations\n",
                        "\tNo missing values allowed in 'A' or 'W'\n",
                        "\t'A' must be binary (0,1)\n",
                        "\t'Z' must be binary (0,1)\n",
                        "\tIntermediate variable (Z) not allowed when there is no experimentation in A",
                        "\tInvalid regression formula for 'Qform'",
                        "\tInvalid regression formula for 'gform'",
                        "\tInvalid regression formula for 'g.Zform'",
                        "\tInvalid regression formula for 'g.Deltaform'",
                        "\tInvalid term name in regression formula for 'Qform'",
                        "\tInvalid term name in regression formula for 'gform'",
                        "\tInvalid term name in regression formula for 'g.Zform'",
                        "\tInvalid term name in regression formula for 'g.Deltaform'"
  )
  if(!all(ok)){
    warning("\n", warning_messages[!ok], immediate. = TRUE)
  }
  return(all(ok))
}


#---------- function .initStage1 ---------------
# Bound Y, map to Ystar if applicable, and
# set boundson on Q and enforce on user-specified values
# returns
#   Ystar - outcome values (between [0,1] if maptoYstar=TRUE)
#   Q - matrix of user-specified values
#   Qbounds - bounds on predicted values for Q (10% wider at each end then
# 			observed range of Y
#			(-Inf,+Inf) is default for linear regression
#   ab - bounding levels used to transform Y to Ystar
#-----------------------------------------------
.initStage1 <- function(Y,A, Q, Q.Z1=NULL, Delta, Qbounds, alpha, maptoYstar, family){
  if(family=="binomial") {Qbounds <- c(0,1)}
  if(is.null(Qbounds)) {
    if(maptoYstar){
      Qbounds <- range(Y[Delta==1])
      Qbounds <- Qbounds + .1*c(-abs(Qbounds[1]),abs(Qbounds[2]))
    } else {
      Qbounds <- c(-Inf, Inf)
    }
  }
  if(!is.null(Q)){
    QAW <- (1-A)*Q[,1] + A*Q[,2]
    Q <- cbind(QAW, Q0W=Q[,1], Q1W=Q[,2])
  }
  if(!is.null(Q.Z1)){
    Q <- cbind(Q, Q0W.Z1=Q.Z1[,1], Q1W.Z1=Q.Z1[,2])
  }
  ab <- c(0,1)
  Ystar <- Y
  if(maptoYstar){
    Ystar <- .bound(Y, Qbounds)
    if(!is.null(Q)){
      Q <- .bound(Q, Qbounds)
    }
    if(0 >= alpha | 1 <= alpha){
      alpha <- .995
      warning(paste("\n\talpha must be between 0 and 1, alpha reset to",alpha,"\n"),
              immediate. = TRUE)
    }
    ab <- range(Ystar, na.rm=TRUE)
    Ystar[is.na(Ystar)] <- 0
    Ystar <- (Ystar-ab[1])/diff(ab)
    if(!is.null(Q)){Q <- (Q-ab[1])/diff(ab)}
    Qbounds <- c(alpha, 1-alpha)
  }
  return(list(Ystar=Ystar, Q=Q, Qbounds=Qbounds, ab=ab))
}

#---------- function .bound ---------------
# set outliers to min/max allowable values
# assumes x contains only numerical data
#-----------------------------------------
.bound <- function(x, bounds){
  x[x>max(bounds)] <- max(bounds)
  x[x<min(bounds)] <- min(bounds)
  return(x)
}


#-----------estimateQ----------------
# purpose: estimate Q=E(Y |Z, A,W) data-adaptively,
# unless super learner not available, or user specifies
# initial values or a regression formula
# arguments:
# 	Y - outcome
# 	Z - intermediate variable between A and Y (default= 0 when no int. var.)
#	A - treatment indicator (1=treatment, 0=control)
# 	W - baseline covariates
#	Delta - missingness indicator
#	Q - optional externally estimated values for Q
#	Qbounds - bounds for predicted values
#  	Qform - optional regression formula to use for glm if
#	        non-data adaptive estimation specified
# 	maptoYstar - if TRUE, using logistic fluctuation for bounded, continuous outcomes
# 		estimation inital Q on linear scale, bounded by (0,1),and return on logit scale
#		(will work if family=poisson)
#	SL.library - library of prediction algorithms for Super Learner
#   cvQinit - flag, if TRUE, cross-validate SL.
# 	family - regression family
#	id - subject identifier
# V - number of cross-validation folds
# returns matrix of linear predictors for Q(A,W), Q(0,W), Q(1,W),
#   (for controlled direct effect, 2 additional columns: Q(Z=1,A=0,W), Q(Z=1,A=1,W))
#		family for stage 2 targeting
#		coef, NA, unless Q is estimated using a parametric model
# 		type, estimation method for Q
#----------------------------------------
estimateQ <- function (Y,Z,A,W, Delta, Q, Qbounds, Qform, maptoYstar,
                       SL.library, cvQinit, family, id, V, verbose) {
  SL.version <- 2
  Qfamily <- family
  m <- NULL
  coef <- NA
  CDE <- length(unique(Z)) > 1
  type <- "user-supplied values"
  if(is.null(Q)){
    if(verbose) { cat("\tEstimating initial regression of Y on A and W\n")}
    Q <- matrix(NA, nrow=length(Y), ncol = 5)
    colnames(Q)<- c("QAW", "Q0W", "Q1W", "Q0W.Z1", "Q1W.Z1")
    if(!(is.null(Qform))){
      if(identical(as.character(stats::as.formula(Qform)), c("~","Y", "."))){
        if(CDE){
          Qform <- paste("Y~Z+A+", paste(colnames(W), collapse="+"))
        } else {
          Qform <- paste("Y~A+", paste(colnames(W), collapse="+"))
        }
      }
      m <- suppressWarnings(stats::glm(Qform, data=data.frame(Y,Z,A,W, Delta), family=family, subset=Delta==1))
      Q[,"QAW"] <- stats::predict(m, newdata=data.frame(Y,Z,A,W), type="response")
      Q[,"Q0W"] <- stats::predict(m, newdata=data.frame(Y,Z=0,A=0,W), type="response")
      Q[,"Q1W"] <- stats::predict(m, newdata=data.frame(Y,Z=0,A=1,W), type="response")
      Q[,"Q0W.Z1"] <- stats::predict(m, newdata=data.frame(Y,Z=1,A=0,W), type="response")
      Q[,"Q1W.Z1"] <- stats::predict(m, newdata=data.frame(Y,Z=1,A=1,W), type="response")
      coef <- stats::coef(m)
      type="glm, user-supplied model"
    } else {
      if(cvQinit){
        m <- try(estQcvSL(Y,X=cbind(Z,A,W),SL.library, family=family,
                          Delta=Delta, Qbounds=Qbounds,id=id, V_SL = V, verbose=verbose))
        if(!(identical(class(m), "try-error"))){
          type <- "cross-validated SL"
          Qinit <- m
          Q <- Qinit$Q
        }
      } else {
        if(verbose) {cat("\t using SuperLearner\n")}
        n <- length(Y)
        X <- data.frame(Z,A,W)
        X00 <- data.frame(Z=0,A=0, W)
        X01 <- data.frame(Z=0,A=1, W)
        newX <- rbind(X, X00, X01)
        if(CDE) {
          X10 <- data.frame(Z=1,A=0, W)
          X11 <- data.frame(Z=1,A=1, W)
          newX <- rbind(newX, X10, X11)
        }
        if(utils::packageDescription("SuperLearner")$Version < SL.version){
          arglist <- list(Y=Y[Delta==1],X=X[Delta==1,], newX=newX, SL.library=SL.library,
                          V=V, family=family, save.fit.library=FALSE, id=id[Delta==1])
        } else {
          arglist <- list(Y=Y[Delta==1],X=X[Delta==1,], newX=newX, SL.library=SL.library,
                          cvControl=list(V=V), family=family, control = list(saveFitLibrary=FALSE), id=id[Delta==1])
        }
        suppressWarnings(
          m<- try(do.call(SuperLearner::SuperLearner, arglist))
        )
        if(identical(class(m),"SuperLearner")){
          Q[,"QAW"] <- m$SL.predict[1:n]
          Q[,"Q0W"] <- m$SL.predict[(n+1):(2*n)]
          Q[,"Q1W"] <- m$SL.predict[(2*n+1):(3*n)]
          if(CDE){
            Q[,"Q0W.Z1"] <- m$SL.predict[(3*n+1):(4*n)]
            Q[,"Q1W.Z1"] <- m$SL.predict[(4*n+1):(5*n)]
          }
          type <- "SuperLearner"
        } else {
          stop("Super Learner failed when estimating Q. Exiting program\n")
        }
      } }
  }
  if(is.na(Q[1,1]) | identical(class(m), "try-error")){
    if(verbose) {cat("\t Running main terms regression for 'Q' using glm\n")}
    Qform <- paste("Y~Z+A+", paste(colnames(W), collapse="+"))
    m <- stats::glm(Qform, data=data.frame(Y,Z,A,W, Delta), family=family, subset=Delta==1)
    Q[,"QAW"] <- stats::predict(m, newdata=data.frame(Y,Z,A,W), type="response")
    Q[,"Q1W"] <- stats::predict(m, newdata=data.frame(Y,Z=0,A=1,W), type="response")
    Q[,"Q0W"] <- stats::predict(m, newdata=data.frame(Y,Z=0,A=0,W), type="response")
    Q[,"Q0W.Z1"] <- stats::predict(m, newdata=data.frame(Y,Z=1,A=0,W), type="response")
    Q[,"Q1W.Z1"] <- stats::predict(m, newdata=data.frame(Y,Z=1,A=1,W), type="response")

    coef <- stats::coef(m)
    type="glm, main terms model"
  }
  Q <- .bound(Q, Qbounds)
  if(maptoYstar | identical(Qfamily,"binomial") | identical(Qfamily, stats::binomial)){
    Q <- stats::qlogis(Q)
    Qfamily <- "binomial"
  } else if (identical(Qfamily, "poisson") | identical(Qfamily, stats::poisson)) {
    Q <- log(Q)
    Qfamily <- "poisson"
  }
  if(!CDE){
    Q <- Q[,1:3]
  }
  if(cvQinit){
    Qinit$Q <- Q
  } else {
    Qinit <- list(Q=Q, family=Qfamily, coef=coef, type=type)
    if(type=="SuperLearner"){
      Qinit$SL.library=SL.library
      Qinit$coef=m$coef
    }
  }
  return(Qinit)
}

#-----------estimateG----------------
# Estimate factors of g
# 		P(A=1|W), P(Z=1|A,W), P(Delta=1|Z,A,W)
# d - dataframe (A,W), (Z,A,W), or (Delta,Z,A,W)
# g1W - optional vector/matrix  of externally estimated values
# gform - optionalformula to use for glm
# SL.library - algorithms to use for super learner estimation
#   id - subject identifier
#  V - number of cross-validation folds
# verbose - flag, whether or not to print messages
# message - printed when verbose=TRUE
# outcome - "A" for treatment, "Z" for intermediate variable,
#           "D" for Delta (missingness)
# newdata - optional values to predict on (needed by tmleMSM function)
# d = [A,W] for treatment
# d = [Z,A,W] for intermediate
# d = [Delta, Z,A,W for missingness]
#----------------------------------------
estimateG <- function (d,g1W, gform, SL.library, id, V, verbose, message, outcome="A", newdata=d)  {
  SL.version <- 2
  SL.ok <- FALSE
  m <- NULL
  coef <- NA
  type <- NULL
  if (is.null(g1W)){
    if(verbose){cat("\tEstimating", message, "\n")}
    if (length(unique(d[,1]))==1) {
      g1W <- rep(1,nrow(d))
      type <- paste("No", strsplit(message, " ")[[1]][1])
      if(outcome=="Z"){
        g1W <- cbind(A0=g1W, A1=g1W)
      } else if (outcome=="D"){
        g1W <- cbind(Z0A0=g1W, Z0A1=g1W, Z1A0=g1W, Z1A1=g1W)
      }
    } else {
      if (is.null(gform)){
        SL.ok <- TRUE
        old.SL <- utils::packageDescription("SuperLearner")$Version < SL.version
        if(old.SL){
          arglist <- list(Y=d[,1], X=d[,-1, drop=FALSE], newX=newdata[,-1, drop=FALSE], family="binomial", SL.library=SL.library, V=V, id=id)
        } else {
          arglist <- list(Y=d[,1], X=d[,-1, drop=FALSE], newX=newdata[,-1, drop=FALSE], family="binomial", SL.library=SL.library, cvControl=list(V=V), id=id)
        }
        suppressWarnings(
          m <- try(do.call(SuperLearner::SuperLearner,arglist))
        )
        if(identical(class(m),"SuperLearner")) {
          g1W <- as.vector(m$SL.predict)
        } else {
          SL.ok <- FALSE
          cat("Error estimating g using SuperLearner. Defaulting to glm\n")
        }
        if (!SL.ok){
          if(verbose){cat("\tRunning main terms regression for 'g' using glm\n")}
          form <- paste(paste(colnames(d)[1],"~1"), paste(colnames(d)[-1], collapse = "+"), sep="+")
          m <- stats::glm(form, data=d, family="binomial")
          g1W <- stats::predict(m, newdata=newdata, type="response")
          coef <- stats::coef(m)
        }
      } else {
        form <- try(stats::as.formula(gform))
        if (inherits(form, "formula")) {
          m <- try(stats::glm(form,  data=d, family="binomial"))
          if (inherits(m, "try-error")){
            if(verbose){cat("\tInvalid formula supplied. Running glm using main terms\n")}
            form <- paste(colnames(d)[1],"~1 + ", paste(colnames(d)[-1], collapse = "+"), sep="")
            m <- stats::glm(form, data=d, family="binomial")
          } else {
            type <- "user-supplied regression formula"
          }
        } else {
          if(verbose){cat("\tRunning main terms regression for 'g' using glm\n")}
          form <- paste(colnames(d)[1],"~1", paste(colnames(d)[-1], collapse = "+"), sep="+")
          m <- stats::glm(form, data=d, family="binomial")
        }
        g1W <- stats::predict(m, newdata=newdata, type="response")
        coef <- stats::coef(m)
      }
      # Get counterfactual predicted values
      if(outcome=="Z"){
        if(identical(class(m),"SuperLearner")){
          g1W <- cbind(stats::predict(m, newdata=data.frame(A=0, newdata[,-(1:2), drop=FALSE]), type="response",
                                      X=d[,-1, drop=FALSE], Y=d[,1])[[1]],
                       stats::predict(m, newdata=data.frame(A=1, newdata[,-(1:2), drop=FALSE]), type="response",
                                      X=d[,-1, drop=FALSE], Y=newdata[,1])[[1]])
        } else {
          g1W <- cbind(stats::predict(m, newdata=data.frame(A=0, newdata[,-(1:2), drop=FALSE]), type="response"),
                       stats::predict(m, newdata=data.frame(A=1, newdata[,-(1:2), drop=FALSE]), type="response"))
        }
        colnames(g1W) <- c("A0", "A1")

      } else if (outcome=="D"){
        if(identical(class(m),"SuperLearner")){
          g1W <- cbind(stats::predict(m, newdata=data.frame(Z=0, A=0, newdata[,-(1:3), drop=FALSE]), type="response",
                                      X=d[,-1,drop=FALSE], Y=d[,1])[[1]],
                       stats::predict(m, newdata=data.frame(Z=0, A=1, newdata[,-(1:3), drop=FALSE]), type="response",
                                      X=d[,-1, drop=FALSE], Y=d[,1])[[1]],
                       stats::predict(m, newdata=data.frame(Z=1, A=0, newdata[,-(1:3), drop=FALSE]), type="response",
                                      X=d[,-1, drop=FALSE], Y=d[,1])[[1]],
                       stats::predict(m, newdata=data.frame(Z=1, A=1, newdata[,-(1:3), drop=FALSE]), type="response",
                                      X=d[,-1, drop=FALSE], Y=d[,1])[[1]])
        } else{
          g1W <- cbind(stats::predict(m, newdata=data.frame(Z=0, A=0, newdata[,-(1:3), drop=FALSE]), type="response"),
                       stats::predict(m, newdata=data.frame(Z=0, A=1, newdata[,-(1:3), drop=FALSE]), type="response"),
                       stats::predict(m, newdata=data.frame(Z=1, A=0, newdata[,-(1:3), drop=FALSE]), type="response"),	     	     		  stats::predict(m, newdata=data.frame(Z=1, A=1, newdata[,-(1:3), drop=FALSE]), type="response"))
        }
        colnames(g1W) <- c("Z0A0", "Z0A1", "Z1A0", "Z1A1")
      }
    }
  } else {
    type <- "user-supplied values"
    if(outcome=="Z") {
      colnames(g1W) <- c("A0", "A1")
    } else if (outcome=="D"){
      colnames(g1W) <- c("Z0A0", "Z0A1", "Z1A0", "Z1A1")[1:ncol(g1W)]
    }
  }
  if(is.null(type)){ type <- class(m)[1]}
  returnVal <- list(g1W=g1W, coef=coef, type=type)
  if(type=="SuperLearner"){
    returnVal$SL.library=SL.library
    returnVal$coef=m$coef
  }
  return(returnVal)
}

# Computes the variance/covariance matrix for all params in the IC.
calcSigma <- function(hAV, gAVW, Y, Q, mAV, covar.MSM, covar.MSMA0, covar.MSMA1, I.V, Delta, ub, id, family, return_IC = FALSE){
  D <- I.V * Delta * ( .bound(hAV[,"hAV"]/gAVW, c(0, ub)) * (Y-Q[,"QAW"]) * covar.MSM
                       + hAV[,"h0V"] * (Q[,"Q0W"] - mAV[,"m0V"]) * covar.MSMA0
                       + hAV[,"h1V"] * (Q[,"Q1W"] - mAV[,"m1V"]) * covar.MSMA1)
  if(!any(is.na(D))){
    # Construct normalizing matrix M
    nterms <- ncol(covar.MSMA0)
    f <- function(x){x[1:nterms] %*% t(x[(nterms+1):(2*nterms)])}
    if(family == "binomial"){
      derivFactor <- cbind(mAV[,"m0V"] * (1-mAV[,"m0V"]), mAV[,"m1V"] * (1-mAV[,"m1V"]))
    } else {
      derivFactor <- matrix(1, nrow=nrow(mAV), ncol = 2)
    }
    deriv.term2 <-  apply(cbind(-hAV[,"h0V"]* I.V * derivFactor[,1] * covar.MSMA0, covar.MSMA0), 1, f)
    deriv.term3 <-  apply(cbind(-hAV[,"h1V"]* I.V * derivFactor[,2] * covar.MSMA1, covar.MSMA1), 1, f)
    ddpsi.D <- as.matrix(deriv.term2 + deriv.term3)
    M <- -matrix(rowMeans(ddpsi.D), nrow=nterms)

    Minv <- try(solve(M))
    if(identical(class(Minv), "try-error")){
      warning("Inference unavailable: normalizing matrix not invertible\n")
      IC    <- NA
      sigma <- NA
    } else {
      Dstar <- t(Minv %*% t(D))
      if(length(unique(id)) < length(Y)){
        Dstar <- matrix(unlist(by(Dstar, id, colMeans, simplify=TRUE)), byrow=TRUE, nrow=length(unique(id)))
      }
      ## Place One
      IC <- Dstar
      sigma <- stats::var(Dstar)
    }
    rownames(sigma) <- colnames(sigma) <- colnames(covar.MSMA0)
  } else {
    D <- .bound(hAV[,"hAV"]/gAVW, c(0, ub)) * (Y-Q[,"QAW"]) * covar.MSM
    term2 <-  hAV[,"h0V"]* (Q[,"Q0W"] - mAV[,"m0V"]) * covar.MSMA0
    term3 <- hAV[,"h1V"]* (Q[,"Q1W"] - mAV[,"m1V"]) * covar.MSMA1
    if (!any(is.na(term2))){D <- D + term2}
    if (!any(is.na(term3))){D <- D + term3}
    ## Place Two
    IC <- I.V * Delta * D
    sigma <- stats::var(IC)
  }
  if(return_IC){
    return(list("sigma" = sigma, "IC" = IC))
  }else{
    return(sigma)
  }
}

my_tmleMSM <- function(Y,A,W,V,T=rep(1,length(Y)), Delta=rep(1, length(Y)), MSM, v=NULL,
                       Q=NULL, Qform=NULL,
                       Qbounds=c(-Inf, Inf),
                       Q.SL.library=c("SL.glm"),
                       cvQinit = FALSE,
                       hAV=NULL, hAVform=NULL,
                       g1W = NULL, gform=NULL,
                       pDelta1=NULL, g.Deltaform=NULL,
                       g.SL.library=c("SL.glm"), ub = 1/0.025,
                       family="gaussian", fluctuation="logistic", alpha  = 0.995,
                       id=1:length(Y), V_SL=5, inference=TRUE, verbose=FALSE, ret_IC = TRUE) {
  Y[is.na(Y)] <- 0
  n <- length(Y)
  n.id <- length(unique(id))
  if(is.null(v)){
    I.V <- rep(1, length(Y))
  } else {
    I.V <- as.numeric(V==v)
  }
  V <- as.matrix(V)
  colnames(V) <- .setColnames(colnames(V), NCOL(V), "V")

  if (sum(!sapply(W, is.numeric)) > 0) {
    stop("Currently, only numeric variables are allowed.  Please convert any character or factor variables to numeric.")
  }

  W <- as.matrix(W)
  colnames(W) <- .setColnames(colnames(W), NCOL(W), "W")

  if(identical (family, stats::binomial)){
    family <- "binomial"
  } else if (identical(family, stats::gaussian)){
    family <- "gaussian"
  }
  if(!.verifyArgs(Y,Z=NULL,A,cbind(V,W,T),Delta, Qform, gform, hAVform, g.Deltaform)){
    stop()
  }
  maptoYstar <- fluctuation=="logistic" | family=="binomial"

  #---- Stage 1 -----
  stage1 <- .initStage1(Y=Y,A=A, Q=Q, Delta=Delta, Qbounds=Qbounds, alpha=alpha, maptoYstar=maptoYstar, family=family)
  Qinit <- suppressWarnings(estimateQ(Y=stage1$Ystar,Z=rep(1, length(Y)), A=A,
                                      W=cbind(W,V,T), Delta=(I.V==1 & Delta==1),
                                      Q=stage1$Q, Qbounds=stage1$Qbounds, Qform=Qform, maptoYstar = maptoYstar,
                                      SL.library=Q.SL.library, cvQinit=cvQinit, family=family, id=id, V = V_SL, verbose=verbose))

  #---- Stage 2 -----
  if(is.null(hAV)){
    gAV <- suppressWarnings(estimateG(d=data.frame(A,V,T), hAV, hAVform,g.SL.library, id, V=V_SL, verbose,
                                      message="h(A,V)", outcome="A"))
    hAV <- cbind((1-A)*(1-gAV$g1W) + A*gAV$g1W, 1-gAV$g1W, gAV$g1W)
  } else {
    hAV <- cbind((1-A)*(hAV[,1]) + A*hAV[,2], hAV)
    gAV <- NULL
    gAV$g1W <- hAV
    gAV$type <- "User-supplied values"
    gAV$coef <- NA
  }
  colnames(hAV) <- c("hAV", "h0V", "h1V")
  if (is.null(v)){
    g <- suppressWarnings(estimateG(d=data.frame(A,V,W,T), g1W, gform,g.SL.library, id, V=V_SL, verbose,
                                    message="treatment mechanism", outcome="A"))
  } else {
    g <- suppressWarnings(estimateG(d=data.frame(A,V,W,T), g1W, gform,g.SL.library, id, V=V_SL, verbose,
                                    message="treatment mechanism", outcome="A", newdata=data.frame(A,V=v, W,T)))
  }
  g$bound <- c(0,ub)
  if(g$type=="try-error"){
    stop("Error estimating treatment mechanism (hint: only numeric variables are allowed)")
  }
  g.Delta <- estimateG(d=data.frame(Delta, Z=1, A, W,V,T), pDelta1, g.Deltaform,
                       g.SL.library,id=id, V = V_SL, verbose, "missingness mechanism", outcome="D")
  g1VW <- g$g1W * g.Delta$g1W[,"Z0A1"]
  g0VW <- (1-g$g1W) * g.Delta$g1W[,"Z0A0"]
  gAVW <- A*g1VW + (1-A)*g0VW

  MSMformula <- stats::formula(paste("Y~",MSM))

  mfA <- stats::model.frame(MSMformula, data=data.frame(Y,A,V,W,T))  # Added Y Oct 24 ,2011
  mf0 <- stats::model.frame(MSMformula, data=data.frame(Y,A=rep(0,n),V,W,T))
  mf1 <- stats::model.frame(MSMformula, data=data.frame(Y,A=rep(1,n),V,W,T))

  covar.MSM   <- stats::model.matrix(MSMformula, mfA)
  covar.MSMA0 <- stats::model.matrix(MSMformula, mf0)
  covar.MSMA1 <- stats::model.matrix(MSMformula, mf1)

  if(verbose){cat("\tTargeting Q\n")}
  C1 <- I.V * .bound(hAV[,"hAV"]/gAVW, c(0,ub)) * covar.MSM
  suppressWarnings(
    epsilon <- stats::coef(stats::glm(stage1$Ystar ~ -1 + offset(Qinit$Q[,"QAW"]) + C1, subset=(I.V==1 & Delta==1), family=Qinit$family))
  )
  Qstar <- cbind(Qinit$Q[,"QAW"] + C1 %*% epsilon,
                 Qinit$Q[,"Q0W"] + .bound(hAV[,"h0V"]/g0VW, c(0,ub)) * covar.MSMA0 %*% epsilon,
                 Qinit$Q[,"Q1W"] + .bound(hAV[,"h1V"]/g1VW, c(0,ub)) * covar.MSMA1 %*% epsilon)
  if(identical(Qinit$family, "binomial")){
    Qstar <- stats::plogis(Qstar)*diff(stage1$ab)+stage1$ab[1]
    Qinit$Q <- stats::plogis(Qinit$Q)*diff(stage1$ab)+stage1$ab[1]
    Ystar <- stage1$Ystar*diff(stage1$ab)+stage1$ab[1]
  }
  colnames(Qstar) <- c("QAW", "Q0W", "Q1W")

  if(verbose){cat("\tEvaluating MSM parameters\n")}
  d.Qstar <- data.frame(Y=c(Qstar[,"Q0W"], Qstar[,"Q1W"]),
                        rbind(mf0, mf1),
                        wts=c(hAV[,"h0V"], hAV[,"h1V"]))
  suppressWarnings(
    psi.Qstar <- stats::coef(stats::glm(MSMformula, data=d.Qstar,
                                        weights=d.Qstar$wts, family=family))
  )
  d.Qinit <- replace(d.Qstar,1, c(Qinit$Q[,"Q0W"], Qinit$Q[,"Q1W"]))
  suppressWarnings(
    psi.Qinit <- stats::coef(stats::glm(MSMformula, data=d.Qinit, weights=d.Qinit$wts, family=family))
  )

  if(inference){
    if(verbose){cat("\tCalculating variance-covariance matrix\n")}
    if(family=="binomial"){
      mAV <- stats::plogis(cbind(covar.MSMA0 %*% psi.Qstar, covar.MSMA1 %*% psi.Qstar))
    } else {
      mAV <- cbind(covar.MSMA0 %*% psi.Qstar, covar.MSMA1 %*% psi.Qstar)
    }
    colnames(mAV) <- c("m0V", "m1V")
    sigma <- calcSigma(hAV, gAVW, Y, Qstar, mAV, covar.MSM, covar.MSMA0, covar.MSMA1, I.V, Delta, ub, id, family,
                       return_IC = ret_IC)
    if(ret_IC){
      IC <- sigma$IC
      sigma <- sigma$sigma/n.id
    }else{
      sigma <- sigma/n.id
      IC <- NULL
    }
    se <- sqrt(diag(sigma))
    pvalue <- 2*stats::pnorm(-abs(psi.Qstar/se))
    lb <- psi.Qstar -1.96 * se
    ub <- psi.Qstar +1.96 * se
  } else {
    sigma <- se <- lb <- ub <- pvalue <- NULL
  }
  Qinit$Q <- Qinit$Q[,-1]
  if(ret_IC){
    returnVal <- list(psi=psi.Qstar, sigma=sigma,se=se, pvalue=pvalue, lb=lb, ub=ub, epsilon=epsilon,  psi.Qinit=psi.Qinit,  Qstar=Qstar[,-1], Qinit=Qinit, g=g, g.AV=gAV, g.Delta=g.Delta, IC=IC)
  }else{
    returnVal <- list(psi=psi.Qstar, sigma=sigma,se=se, pvalue=pvalue, lb=lb, ub=ub, epsilon=epsilon,  psi.Qinit=psi.Qinit,  Qstar=Qstar[,-1], Qinit=Qinit, g=g, g.AV=gAV, g.Delta=g.Delta)
  }
  class(returnVal) <- "tmleMSM"
  return(returnVal)
}

