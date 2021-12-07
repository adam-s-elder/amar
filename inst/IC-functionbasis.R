#' Functional Null estimate based on the work of Ted Westling
#' and suggestions from Yen-Chi Chen
#'
#' @param obs_dat The observed data.
#' @param what should the IC, estimate, or both be returned.
#' @param control other control parameters.
#'
#' @export

ic.fb.est <- function(obs_dat, what = "both", control = NULL) {
  Y <- obs_dat[, 1, drop = TRUE]
  A <- obs_dat[, 2, drop = TRUE]
  if (!is.null(control$cross.fit.func)) {
    crsfit <- control$cross.fit.func
    V <- 5
  }else{
    crsfit <- FALSE
    V <- NULL
  }
  W <- data.frame(W1 = obs_dat[, 3, drop = TRUE])
  cat("\n Starting function estimation \n")
  tictoc::tic()
  both <- causalNullTest(Y, A, W, p = NA,
                         control = list(cross.fit = crsfit, V = V,
                                        verbose = FALSE,
                                        g.n.bins = 2:5),
                         return_ic_est = TRUE)

  cat("\n Finished with function estimation \n", do.call(tictoc::toc.outmsg,
                                                         tictoc::toc()), " \n")
  tictoc::tic()
  cat("\n Starting Basis conversion \n")
  if (!is.null(control$num_basis)) {
    both <- convert_to_basis(est = both$est, ic = both$ic,
                             gridvals = both$gridpnts,
                             num_basis = control$num_basis)
  } else {
    both <- convert_to_basis(est = both$est, ic = both$ic,
                             gridvals = both$gridpnts)
  }
  cat("\n Done with basis conversion \n", do.call(tictoc::toc.outmsg,
                                                  tictoc::toc()), " \n")

  fin <- list()
  if (what %in% c("both", "est")) {
    fin$est <- both$est
  }
  if (what %in% c("both", "ic")) {
    fin$ic <- both$ic
  }
  return(fin)
}

#' @export

causalNullTest <- function(Y, A, W, p=2,
                           control = list(),
                           return_ic_est = FALSE) {

  call <- match.call(expand.dots = TRUE)
  control <- do.call("causalNullTest.control", control)

  .check.input(Y = Y, A = A, W = W, p = p, control = control)

  require(mvtnorm)
  n <- length(Y)
  a.vals <- sort(unique(A))

  if(control$cross.fit & is.null(control$folds)) {
    control$folds <- sample(rep(1:control$V, length.out = n), replace=FALSE)
  }

  if(is.null(control$mu.hat)) {
    library(SuperLearner)
    if(control$cross.fit) {
      if(control$verbose) cat("Estimating outcome regressions...")
      control$mu.hat <- lapply(1:control$V, function(v) {
        if(control$verbose) cat("fold", v, "...")
        if(length(setdiff(Y, c(0,1))) == 0) {
          fit <- SuperLearner(Y = Y[control$folds != v], X = cbind(A, W)[control$folds != v,],
                              SL.library = control$mu.SL.library,
                              family = 'binomial',
                              method = 'method.NNloglik')
        } else {
          fit <- SuperLearner(Y = Y[control$folds != v], X = cbind(A, W)[control$folds != v,],
                              SL.library = control$mu.SL.library,
                              family = 'gaussian', method = 'method.CC_LS')
        }
        function(a, w) c(predict(fit, newdata = cbind(A = a, w),onlySL = TRUE)$pred)
      })
      if(control$verbose) cat("\n")
    } else {
      if(control$verbose) cat("Estimating outcome regression...")
      if(length(setdiff(Y, c(0,1))) == 0) {
        mu.fit <- SuperLearner(Y = Y, X = cbind(A, W), SL.library = control$mu.SL.library,
                               family = 'binomial', method = 'method.NNloglik')
      } else {
        mu.fit <- SuperLearner(Y = Y, X = cbind(A, W), SL.library = control$mu.SL.library,
                               family = 'gaussian', method = 'method.CC_LS')
      }
      control$mu.hat <- function(a, w) c(predict(mu.fit, newdata = cbind(A = a, w),onlySL = TRUE)$pred)
      if(control$verbose) cat("\n")
    }
  }

  if(is.null(control$g.hat)) {
    if(control$cross.fit) {
      if(control$verbose) cat("Estimating propensities..")
      control$g.hat <- lapply(1:control$V, function(v) {
        if(control$verbose) cat("fold", v)
        fit <- cmdSuperLearner(A = A[control$folds != v], W = W[control$folds != v,,drop=FALSE], newA = A[control$folds == v], newW = W[control$folds == v,,drop=FALSE], control=list(SL.library = control$g.SL.library, n.bins = control$g.n.bins, verbose = control$verbose, saveFitLibrary = FALSE))
        c(fit$SL.densities)
        #function(a, w) c(predict.cmdSuperLearner(fit, newA = a, newW = w))
      })
      if(control$verbose) cat("\n")
    } else {
      if(control$verbose) cat("Estimating propensity...")
      g.fit <- cmdSuperLearner(A = A, W = W,
                               control=list(SL.library = control$g.SL.library,
                                            n.bins = control$g.n.bins,
                                            verbose = control$verbose,
                                            saveFitLibrary = FALSE))
      control$g.hat <- g.fit$SL.densities
      rm(g.fit)
      #control$g.hat <- function(a, w) c(predict.cmdSuperLearner(g.fit, newA = a, newW = w))

      if(control$verbose) cat("\n")
    }
  }

  if(control$verbose) cat("Computing Omega...")
  if(!control$cross.fit) {
    ord <- order(A)
    A <- A[ord]
    Y <- Y[ord]
    W <- W[ord,,drop=FALSE]
    if(inherits(control$g.hat, "function")) {
      g.hats <- control$g.hat(A, W)
      control$g.hat <- NULL
    }
    else g.hats <- control$g.hat
    if(any(g.hats < control$g.trunc)) {
      warning("Truncating g.hats below. Possible positivity issues.")
      g.hats[g.hats < control$g.trunc] <- control$g.trunc
    }
    a.ecdf <- ecdf(A)
    a.weights <- sapply(a.vals, function(a0) mean(A == a0))
    A.a.val <- sapply(A, function(a0) which(a.vals == a0))
    u.vals <- a.ecdf(a.vals)
    mu.hats.a.vals <- sapply(a.vals, function(a0) control$mu.hat(a0, W)) #rows index W, columns index a.vals
    control$mu.hat <- NULL
    mu.hats <- mu.hats.a.vals[,A.a.val]
    theta.a.vals <- colMeans(mu.hats.a.vals)
    theta.A <- theta.a.vals[A.a.val]
    mu.hats.data <- diag(mu.hats)
    partial.mu.means <- t(apply(mu.hats, 1, cumsum)) / n
    gamma.hat <- mean(mu.hats)
    Omega.a.vals <- sapply(a.vals, function(a0) mean(as.numeric(A <= a0) * theta.A)) - gamma.hat * u.vals

    IF.vals <- sapply(a.vals, function(a0) {
      if(any(A <= a0))  mumean.vals <- partial.mu.means[,max(which(A <= a0))]
      else mumean.vals <- 0
      (as.numeric(A <= a0) - a.ecdf(a0)) * ((Y - mu.hats.data) / g.hats + theta.A - gamma.hat) + mumean.vals - partial.mu.means[,n] * a.ecdf(a0) - 2 * Omega.a.vals[which(a.vals == a0)]
    })

    Omega.hat <- colMeans(IF.vals) + Omega.a.vals

    if(control$verbose) cat("\nComputing covariance...\n")

    if (return_ic_est){
      return(list("est" = Omega.hat, "ic" = IF.vals, "gridpnts" = a.vals))
    }else{
      Sigma.hat <- sapply(1:length(a.vals), function(s) sapply(1:length(a.vals), function(t) {
        mean(IF.vals[,s] * IF.vals[,t])
      }))
    }
  }else{
    fold.Omega.hats <- matrix(NA, nrow = control$V, ncol = length(a.vals))
    IF.vals <- vector(length=control$V, mode='list')
    for(j in 1:control$V) {
      if(control$verbose) cat("fold", j, "...")
      Nv <- sum(control$folds == j)
      A.test <- A[control$folds == j]
      Y.test <- Y[control$folds == j]
      W.test <- W[control$folds == j,, drop=FALSE]
      ord <- order(A.test)
      A.test <- A.test[ord]
      Y.test <- Y.test[ord]
      W.test <- W.test[ord,, drop=FALSE]
      if(inherits(control$g.hat[[j]], "function")) {
        g.hats.test <- control$g.hat[[j]](a = A.test, w = W.test)
      }
      else g.hats.test <- control$g.hat[[j]]

      if(any(g.hats.test < control$g.trunc)) {
        warning("Truncating g.hats below. Possible positivity issues.")
        g.hats.test[g.hats.test < control$g.trunc] <- control$g.trunc
      }
      a.ecdf <- ecdf(A.test)
      a.weights <- sapply(a.vals, function(a0) mean(A.test == a0))
      A.a.val <- sapply(A.test, function(a0) which(a.vals == a0))
      u.vals <- a.ecdf(a.vals)
      mu.hats.a.vals <- sapply(a.vals, function(a0) control$mu.hat[[j]](a=a0, w=W.test)) #rows index W, columns index a.vals
      mu.hats <- mu.hats.a.vals[,A.a.val]
      theta.a.vals <- colMeans(mu.hats.a.vals)
      theta.A <- theta.a.vals[A.a.val]
      mu.hats.data <- diag(mu.hats)
      partial.mu.means <- t(apply(mu.hats, 1, cumsum)) / Nv
      gamma.hat <- mean(mu.hats)
      Omega.a.vals <- sapply(a.vals, function(a0) mean(as.numeric(A.test <= a0) * theta.A)) - gamma.hat * u.vals

      IF.vals[[j]] <- sapply(a.vals, function(a0) {
        if(any(A.test <= a0)) mumean.vals <- partial.mu.means[,max(which(A.test <= a0))]
        else mumean.vals <- 0
        (as.numeric(A.test <= a0) - a.ecdf(a0)) * ((Y.test - mu.hats.data) / g.hats.test + theta.A - gamma.hat) + mumean.vals - partial.mu.means[,ncol(partial.mu.means)] * a.ecdf(a0) - 2 * Omega.a.vals[which(a.vals == a0)]
      })
      fold.Omega.hats[j,] <- colMeans(IF.vals[[j]]) + Omega.a.vals
    }
    Omega.hat <- colMeans(fold.Omega.hats)
    if(control$verbose) cat("\nComputing covariance...\n")
    if (return_ic_est) {
      IF.vals <- do.call(rbind, IF.vals)
      return(list("est" = Omega.hat, "ic" = IF.vals))
    }else{
      Sigma.hat <- sapply(1:length(a.vals),
                          function(s) sapply(1:length(a.vals), function(t) {
                            mean(unlist(lapply(IF.vals, function(IF) mean(IF[,s] * IF[,t]))))
                          }))
    }
  }

  if(control$verbose) cat("Simulating paths...\n")

  paths <- rmvnorm(control$n.sim, sigma=Sigma.hat)

  if(control$verbose) cat("Computing statistics...\n")

  a.weights <- sapply(a.vals, function(a) mean(A == a))
  ret <- t(sapply(p, function(pp) {
    stat <- ifelse(pp < Inf, (sum(abs(Omega.hat )^pp * a.weights))^{1/pp}, max(abs(Omega.hat)))

    if(pp < Inf) {
      stats <- (apply(abs(paths)^pp, 1, function(row) sum(row * a.weights)))^{1/pp}
    } else {
      stats <- apply(abs(paths), 1, max)
    }

    p.val <- mean(stats / sqrt(n) > stat)

    q <- quantile(stats, (1 - (1-control$conf.level) / 2))
    ci.ll <- max(stat - q / sqrt(n), 0)
    ci.ul <- stat + q / sqrt(n)

    res <- c(stat, p.val, ci.ll, ci.ul)

    res

  }))
  ret.df <- data.frame(p = p, obs.stat = ret[,1], p.val = ret[,2], ci.ll = ret[,3], ci.ul = ret[,4])
  ret.list <- list(test = ret.df)
  if (control$return.Omega) {
    ret.list <- c(ret.list, list(Omega.hat = data.frame(a = a.vals, Omega.hat),
                                 IF.vals = IF.vals, paths = paths))
  }
  if (control$save.nuis.fits) {
    ret.list <- c(ret.list, mu.hat = control$mu.hat, g.hat = control$g.hat)
    if (control$cross.fit) ret.list <- c(ret.list, folds = control$folds)
  }

  return(ret.list)
}


#' Initialize control parameters for causalNullTest
#'
#' This function initializes the control parameters for use in \code{causalNullTest}. The outcome regression function mu is by default estimated using \code{\link[SuperLearner]{SuperLearner}}, and the propensity is estimated using the conditional mixed density method implemented in \code{\link{cmdSuperLearner}}. Alternatively, the estimation process can be overriden by providing predictions from pre-fit nuisance estimators.
#'
#' @param mu.SL.library Library of candidate learners for the outcome regression to be passed on to \code{\link[SuperLearner]{SuperLearner}}. Ignored if \code{mu.hats} is provided.
#' @param g.SL.library Library of candidate learners for the outcome regression to be passed on to \code{\link{cmdSuperLearner}}.
#' @param g.n.bins Numeric vector of number of bins to use for estimation of the propensity. Passed on to \code{\link{cmdSuperLearner}}.
#' @param cross.fit Logical indicating whether to cross-fit nuisance parameters. Defaults to \code{TRUE}.
#' @param V Positive integer number of folds for cross-fitting. Defaults to 10.
#' @param folds Optional \code{n x 1} numeric vector indicating which fold each observation is in.
#' @param save.nuis.fits Logical indicating whether to save the fitted nuisance objects.
#' @param mu.hat Optional pre-fit outcome regression. If \code{cross.fit} is \code{FALSE}, then a function that takes arguments \code{a} (a vector) and \code{w} (a data.frame) and returns predictions of the outcome regression function. If \code{cross.fit} is \code{TRUE}, then a list of functions of length \code{V} with the fitted outcome regression functions on each of the \code{V} training sets. If provided as a list of functions, then \code{folds} must be provided. If \code{mu.SL.library = NULL}, then \code{mu.hats} must be specified.
#' @param g.hat Optional pre-fit treatment propensities. If \code{cross.fit} is \code{FALSE}, then a function that takes arguments \code{a} (a vector) and \code{w} (a data.frame) and returns predictions of the standardized propensity function. If \code{cross.fit} is \code{TRUE}, then a list of functions of length \code{V} with the fitted standardized propensity functions on each of the \code{V} training sets. If provided as a list of functions, then \code{folds} must be provided. If \code{g.SL.library = NULL}, then \code{g.hats} must be specified.
#' @param g.trunc Value at which to truncate predicted propensities from below. Any propensity values less than \code{g.trunc} will be set to \code{g.trunc}.
#' @param n.sim Number of simulations to use for the limiting Gaussian process in computing approximate quantiles.
#' @param return.Omega Logical indicating whether to return the estimated primitive function Omega.
#' @param conf.level Optional confidence level to use for computing confidence bands for Omega.
#' @param verbose Logical indicating whether to print progress to the command line.
#' @return Named list containing the control options.

causalNullTest.control <- function(mu.SL.library = c("SL.mean", "SL.glm", "SL.gam", "SL.earth"),
                                   g.SL.library = c("SL.mean", "SL.glm", "SL.gam", "SL.earth"),
                                   g.n.bins = 2:(length(unique(A))/50),
                                   cross.fit = TRUE,
                                   V = 10,
                                   folds = NULL,
                                   save.nuis.fits = FALSE,
                                   mu.hat = NULL,
                                   g.hat = NULL,
                                   g.trunc = .001, n.sim = 1e4,
                                   return.Omega = FALSE,
                                   conf.level = .95,
                                   verbose = FALSE) {
  list(mu.SL.library = mu.SL.library, g.SL.library = g.SL.library, g.n.bins = g.n.bins, cross.fit = cross.fit, V = V,
       folds = folds, mu.hat = mu.hat,  g.hat = g.hat, g.trunc = g.trunc, n.sim = n.sim, return.Omega = return.Omega,
       save.nuis.fits = save.nuis.fits, conf.level = conf.level, verbose = verbose)
}

.check.input <- function(Y, A, W, p, control) {
  if(length(Y) != length(A)) stop("Y and A must have the same length")
  if(!is.null(control$mu.hat)) {
    if(control$cross.fit) {
      if(is.null(control$folds)) {
        stop("mu.hat provided and cross.fit=TRUE, but folds not provided.")
      }
      if(length(control$mu.hat) < length(unique(control$folds))) {
        stop("mu.hat provided and cross.fit=TRUE, but mu.hats is not a list with the same length as the number of folds.")
      }
    } else {
      if(!is.null(control$folds)) {
        stop("mu.hat provided and cross.fit=FALSE, but folds were provided.")
      }
    }
  }
  if(!is.null(control$g.hat)) {
    if(control$cross.fit) {
      if(is.null(control$folds)) {
        stop("g.hat provided and cross.fit=TRUE, but folds not provided.")
      }
      if(length(control$g.hat) < length(unique(control$folds))) {
        stop("g.hat provided and cross.fit=TRUE, but g.hats is not a list with the same length as the number of folds.")
      }
    } else {
      if(!is.null(control$folds)) {
        stop("g.hat provided and cross.fit=FALSE, but folds were provided.")
      }
    }
  }
  if(is.null(control$mu.SL.library)) {
    if(is.null(control$mu.hats)) {
      stop("mu.hats must be provided if mu.SL.library is not specified.")
    }
  } else {
    if(is.null(control$mu.hats)) {
      if(control$cross.fit & is.null(control$folds) & is.null(control$V)) {
        stop("cross.fit = TRUE, but number of folds not specified.")
      }
    }
  }
  if(is.null(control$g.SL.library)) {
    if(is.null(control$g.hats)) {
      stop("g.hats must be provided if g.SL.library is not specified.")
    }
  } else {
    if(is.null(control$g.hats)) {
      if(control$cross.fit & is.null(control$folds) & is.null(control$V)) {
        stop("cross.fit = TRUE, but number of folds not specified.")
      }
    }
  }
  if(any(is.na(Y) | is.na(A) | is.na(W))) {
    stop("Missing outcome, treatment, or confounders detected; missing data not allowed.")
  }
}

#' SuperLearner-based estimation of (c)onditional (m)ixed continuous-discrete (d)ensity function
#'
#' This function estimates a standardized conditional density function that may have both continuous and discrete components. Let \code{A} be a univariate exposure and \code{W} be a p-dimensional vector of covariates. Then this function estimates p(a | w) / p(a) at points  of absolute continuity of the marginal distribution of \code{A}, where p(a | w) = (d/da)P(A <= a | W = w) is the conditional density of \code{A} given \code{W = w} evaluated at a and p(a) = (d/da) P(A <= a) is the marginal density of \code{A}, and at discrete points of the marginal distribution of \code{A}, this function estimates P(A = a | W = w)/P(A = a).
#'
#' The basic idea is to first transform A by its empirical CDF to obtain U = F_n(A), because the conditional density or  mass function of F(A) equals the standardized conditional density/mass of A for F(a) = P(A <= a). Then, the support [0,1] of U is discretized into \code{b} sets (which may be singleton sets) using the marginal distribution of U. Within each of these sets, the conditional probability that U falls in the set given \code{W} is estimated using the specified wrapper algorithms from the \code{\link[SuperLearner]{SuperLearner}} package. This procedure is repeated over a set of possible number of bins \code{b}, and optimal weights for all algorithms are found using negative log likelihood loss.
#'
#' @param A \code{n x 1} numeric vector of exposure values.
#' @param W \code{n x p} data.frame of covariate values to condition upon.
#' @param newA \code{m x 1} numeric vector of new exposure values at which to obtain predictions. Defaults to \code{A}.
#' @param newW \code{m x p} data.frame of new covariate values at which to obtain predictions. Defaults to \code{W}.
#' @param control Optional list of control parameters. See \code{\link{cmdSuperLearner.control}} for details.
#' @param cvControl Optional list of control parameters for cross-validation. See \code{\link{cmdSuperLearner.cvControl}} for details.
#' @return \code{cmdSuperLearner} returns a named list with the following elements:
#' \item{fits}{A list of fits for each of the number of bins specified in control$n.bins, as output by \link{cmdSuperLearner.onebin.}}
#' \item{cv.library.densities}{Cross-validated densities from every element of the library.}
#' \item{library.densities}{Densities predicted using the new data.}
#' \item{SL.densities}{Super learner densities predicted on the new data.}
#' \item{coef}{The coefficient of the meta-learner.}
#' \item{library.names}{Names of library algortihms.}
#' \item{a.ecdf}{Empirical CDF of the exposure.}
#' \item{control}{Control elements used in fitting.}
#' \item{cvControl}{Cross-validation controls used in fitting.}
#'
#' @examples
#' # Sample data
#' #n <- 1000
#' #W <- data.frame(W1 = runif(n))
#' #Z <- rbinom(n, size = 1, prob = 1/(1 + exp(2-W$W1)))
#' #A <- (1-Z) * rnorm(n, mean = W$W1, sd = abs(1 + W$W1))
#' #fit <- cmdSuperLearner(A, W, control=list(SL.library = c("SL.mean", "SL.glm", "SL.gam", "SL.earth"), verbose=TRUE, n.bins = c(2:10)))


cmdSuperLearner <- function(A, W, newA = A, newW = W, control = list(), cvControl = list()) {
  n <- nrow(W)

  call <- match.call(expand.dots = TRUE)
  control <- do.call("cmdSuperLearner.control", control)
  cvControl <- do.call("cmdSuperLearner.CV.control", cvControl)

  validRows <- cmdCVFolds(n = n, cvControl = cvControl)

  library(Rsolnp)

  fits <- NULL
  for(b in control$n.bins) {
    if(control$verbose) cat("\nEstimating models with", b, "bins... ")
    fits[[paste0('dens.fit.', b, 'bins')]] <- cmdSuperLearner.onebin(A, W, newA=newA, newW = newW, b=b, SL.library = control$SL.library, verbose = control$verbose, validRows = validRows, saveFitLibrary = control$saveFitLibrary)
  }

  algs.per.bin <- ncol(fits[[1]]$cv.library.densities)
  n.algs <- length(control$n.bins) * algs.per.bin
  cv.library.densities <- matrix(NA, nrow=n, ncol=n.algs)
  library.densities <- matrix(NA, nrow=length(newA), ncol=n.algs)
  library.names <- NULL
  start.col <- 1
  for(b in control$n.bins) {
    end.col <- start.col + algs.per.bin - 1
    cv.library.densities[,start.col:end.col] <- fits[[paste0('dens.fit.', b, 'bins')]]$cv.library.densities
    library.densities[,start.col:end.col] <- fits[[paste0('dens.fit.', b, 'bins')]]$library.densities
    library.names <- c(library.names, fits[[paste0('dens.fit.', b, 'bins')]]$alg.names)
    start.col <- end.col + 1
  }

  if(control$verbose) cat("\nOptimizing model weights...\n")

  # Remove algs with errors in cv predictions
  errors.in.library <- apply(cv.library.densities, 2, function(col) any(is.na(col)))
  if(any(errors.in.library)) warning(paste0("Errors in the following candidate algorithms: ", library.names[which(errors.in.library)]))
  n.include <- sum(!errors.in.library)

  # Do SL log-likelihood optimization
  cv_risk <- function(beta) -mean(log(cv.library.densities[,!errors.in.library] %*% beta))
  capture.output(solnp_solution <- solnp(rep(1/n.include, n.include), cv_risk, eqfun=sum, eqB=1, ineqfun=function(beta) beta, ineqLB=rep(0,n.include), ineqUB=rep(1, n.include)))
  coef <- rep(0, n.algs)
  coef[!errors.in.library] <- solnp_solution$pars
  if(control$verbose) {
    cat("Top five learners by weight: \n")
    for(j in 1:5) {
      cat(library.names[order(coef, decreasing = TRUE)[j]], " (weight ", sort(coef, decreasing = TRUE)[j], ")\n", sep='')
    }
  }
  SL.density <- c(library.densities[,!errors.in.library,drop=FALSE] %*% solnp_solution$pars)

  return(list(fits = fits, cv.library.densities = cv.library.densities, library.densities = library.densities, SL.densities = SL.density, coef = coef, library.names = library.names, a.ecdf = ecdf(A), control=control, cvControl = cvControl))


}


#' Control parameters for the conditional mixed density Super Learner
#'
#' This function initiates control parameters for the \code{\link{cmdSuperLearner}} function.
#'
#' @param n.bins Vector of integers >= 2 indicating the number of bins to use for discretization. Defaults to \code{2:floor(n/100)}.
#' @param SL.library Library to use for bin-specific SuperLearners. Defaults to \code{c("SL.mean", "SL.glm", "SL.gam", "SL.earth")}.
#' @param saveFitLibrary Logical indicating whether to save the fit library on the full data. Defaults to \code{TRUE}. If \code{FALSE}, cannot obtain predicted values on new data later.
#' @param verbose Logical indicating whether to report progress of the estimation process.
#' @return Returns a named list with control parameters.

cmdSuperLearner.control <- function (n.bins = 2:floor(length(unique(A))/50), SL.library = c("SL.mean", "SL.glm", "SL.gam", "SL.earth"), saveFitLibrary = TRUE, verbose = FALSE) {
  list(n.bins = n.bins, SL.library = SL.library, saveFitLibrary = saveFitLibrary, verbose = verbose)
}

#' Control parameters for the cross validation steps in conditional mixed density Super Learner
#'
#' This function initiates control parameters for the cross-validation in \code{\link{cmdSuperLearner}} function.
#'
#' @param V Number of cross-validation folds. Defaults to 10.
#' @param shuffle Logical indicating whether to shuffle the indices, or to simply assign sequentially. Defaults to \code{TRUE}. Should almost always be set to \code{TRUE} unless it is explicitly desired to assign sequentially.
#' @param validRows Optional custom list of indices for validation folds.
#' @return Returns a list of length \code{V} with validation indices for each of the folds.

cmdSuperLearner.CV.control <- function (V = 10L, shuffle = TRUE, validRows = NULL) {
  V <- as.integer(V)
  if (!is.null(validRows)) {
    if (!is.list(validRows)) {
      stop("validRows must be a list of length V containing the row numbers for the corresponding validation set")
    }
    if (!identical(V, length(validRows))) {
      stop("V and length(validRows) must be identical")
    }
  }
  list(V = V, shuffle = shuffle, validRows = validRows)
}

#' Create cross-validation folds
#'
#' This function generates cross-validation folds.
#'
#' @param n Number of observations.
#' @param cvControl Named list generated by \link{cmdSuperLearner.CV.control}.

cmdCVFolds <- function (n, cvControl) {
  if (!is.null(cvControl$validRows)) return(cvControl$validRows)
  stratifyCV <- cvControl$stratifyCV
  shuffle <- cvControl$shuffle
  V <- cvControl$V
  if (shuffle) {
    validRows <- split(sample(1:n), rep(1:V, length = n))
  }
  else {
    validRows <- split(1:n, rep(1:V, length = n))
  }

  return(validRows)
}


.find.bin <- function(x, bins) {
  mat <- t(sapply(x-1e-10, function(x0) {
    unlist(lapply(bins, function(bin) {
      interval_contains_element(bin, x0)
    }))
  }))
  if(any(rowSums(mat) > 1)) stop("Overlapping bins")
  if(any(rowSums(mat) == 0)) stop("Element outside all bins")

  apply(mat, 1, function(row) which(row))
}

#' cmdSuperLearner for a specific number of bins
#'
#' This function estimates the conditional mixed density using a given number of bins \code{b}.
#'
#' @param A \code{n x 1} numeric vector of exposure values.
#' @param W \code{n x p} data.frame of covariate values to condition upon.
#' @param newA \code{m x 1} numeric vector of new exposure values at which to obtain predictions. Defaults to \code{A}.
#' @param newW \code{m x p} data.frame of new covariate values at which to obtain predictions. Defaults to \code{W}.
#' @param b Integer number of bins >= 2.
#' @param SL.library Library to use for bin-specific probabilities.
#' @param verbose Logical indicating whether to print progress reports to the command line.
#' @param validRows List of rows in each CV fold.
#' @return Returns a named list with the following elements:
#' \item{bins}{List of length \code{b} containing the sets used for each bin.}
#' \item{bin.fits}{List of length \code{b} containing the estimated SuperLearner objects for each bin.}
#' \item{a.ecdf}{Empirical CDF of the exposure.}
#' \item{SL.bin.probs}{SuperLearner conditional probabilities of being in each bin for the new data.}
#' \item{SL.densities}{SuperLearner conditional standardized mixed density correspondint to each bin for the new data.}
#' \item{cv.library.densities}{Cross-validated library conditional standardized mixed density corresponding to each bin.}
#' \item{library.densities}{Library conditional standardized mixed density corresponding to each bin fit on the new data.}
#' \item{alg.names}{Algorithm names.}
#' @examples
#' # Define parameters
#' #n <- 300
#' #W <- data.frame(matrix(rnorm(3 * n), ncol = 3))
#' #Z <- rbinom(n, size = 1, prob = 1/(1 + exp(2-W[,1] + W[,2])))
#' #A <- (1-Z) * rnorm(n, mean = W[,2] - W[,3], sd = abs(1 + W[,1]))
#' #validRows <- cmdCVFolds(n = n, cvControl = list(V = 10, shuffle=TRUE, validRows = NULL))
#' #bin.fit <- cmdSuperLearner.onebin(A, W, b = 2, SL.library = c("SL.mean", "SL.glm", "SL.gam", "SL.earth"), verbose=TRUE, validRows = validRows)

cmdSuperLearner.onebin <- function(A, W, newA=A, newW=W, b, SL.library, verbose, validRows, saveFitLibrary) {
  n.folds <- length(validRows)
  a.ecdf <- ecdf(A)
  U <- a.ecdf(A)
  n <- nrow(W)
  m <- nrow(newW)
  W <- as.data.frame(W)
  newW <- as.data.frame(newW)
  U <- as.numeric(U)
  library(Rsolnp)
  library(SuperLearner)
  library(sets)
  tab <- table(U)
  un.U <- as.numeric(names(tab))
  un.U.frac <- as.numeric(tab) / length(U)
  if(b <= 1) stop("Number of bins must be > 1")
  if(length(un.U) < b) stop("Number of bins must not be larger than number of unique values of U.")
  if(length(un.U) == b) {
    mass.pts <- un.U
    bins <- data.frame(bin = 1:b, lower = un.U, upper = un.U, bin.length = 0, mass.pt = TRUE)
  }
  if(length(un.U) > b) {
    if(any(un.U.frac >= 1/b)) {
      mass.pts <- un.U[un.U.frac >= 1/b]
      n.mass.pts <- length(mass.pts)
      mass.pt.lowers <- round(sapply(mass.pts, function(x) max(c(U[x - U > 1/(10*n)], 0))), 7)
      mass.intervals <- lapply(1:n.mass.pts, function(j) {
        interval(mass.pt.lowers[j], mass.pts[j], bounds="(]")
      })
      cont.intervals <- data.frame(lower=c(0,mass.pts), upper=c(mass.pt.lowers, 1))
      cont.intervals$length <- cont.intervals$upper - cont.intervals$lower
      cont.intervals <- subset(cont.intervals, length > 0)
    }
    else {
      mass.pts <- NULL
      n.mass.pts <- 0
      mass.intervals <- NULL
      cont.intervals <- data.frame(lower=0, upper=1, length=1)
    }

    n.cont.bins <- b - n.mass.pts
    if(n.cont.bins > 0) {
      delta <- sum(cont.intervals$length) / n.cont.bins
      delta <- round(delta, digits=ceiling(log10(n)) + 2)
      cont.bin.endpts <- matrix(NA, nrow=n.cont.bins, ncol=2)

      for(j in 1:n.cont.bins) {
        if(j == 1) start <- cont.intervals$lower[1]
        else start <- end
        start.interval <- max(which(cont.intervals$lower <= start + 1e-6 & start <= cont.intervals$upper + 1e-6))
        if(start == cont.intervals$upper[start.interval]) {
          start.interval <- start.interval + 1
          start <- cont.intervals$lower[start.interval]
        }
        end <- start + delta
        end.interval <- start.interval
        if(!(all.equal(end, cont.intervals$upper[end.interval]) == TRUE) && end > cont.intervals$upper[end.interval]) {
          length.used <- cont.intervals$upper[end.interval] - start
          length.left <- delta - length.used
          end.interval <- end.interval + 1
          end <- cont.intervals$lower[end.interval] + length.left
        }
        while(!(all.equal(end, cont.intervals$upper[end.interval]) == TRUE) && end > cont.intervals$upper[end.interval]) {
          length.used <- length.used + cont.intervals$upper[end.interval] - cont.intervals$lower[end.interval]
          length.left <- delta - length.used
          end.interval <- end.interval + 1
          end <- cont.intervals$lower[end.interval] + length.left
        }
        end <- round(end, digits=ceiling(log10(n)) + 3)
        if(j == n.cont.bins) end <- cont.intervals$upper[nrow(cont.intervals)]
        cont.bin.endpts[j,] <- c(start, end)
      }

      cont.intervals <- lapply(1:n.cont.bins, function(j) {
        if(j == 1) int <- interval(cont.bin.endpts[j, 1], cont.bin.endpts[j, 2], bounds="[]")
        else int <- interval(cont.bin.endpts[j, 1], cont.bin.endpts[j, 2], bounds="(]")
        if(n.mass.pts > 0) {
          for(k in 1:n.mass.pts) {
            int <- interval_complement(mass.intervals[[k]], int)
          }
        }
        return(int)
      })
    } else {
      cont.intervals <- list()
    }

    bins <- c(mass.intervals, cont.intervals)
  }

  bin.sizes <- unlist(lapply(bins, interval_measure))

  disc.U <- .find.bin(U, bins)

  U.new <- a.ecdf(newA)
  disc.U.new <- .find.bin(U.new, bins)

  bin.fracs <- sapply(1:b, function(j) mean(disc.U == j))

  bin.fits <- NULL
  bin.probs <- matrix(NA, nrow=n, ncol=b)
  for(bin in 1:b) {
    if(verbose) cat("bin", bin, "... ")
    capture.output(bin.fit <- try(SuperLearner(Y=as.numeric(disc.U==bin), X=W, newX=newW, family='binomial', SL.library = SL.library, method='method.NNloglik', control = list(saveFitLibrary=saveFitLibrary), cvControl = list(V=n.folds, validRows=validRows)), silent=TRUE))
    if(class(bin.fit) == "try-error") {
      capture.output(bin.fit <- try(SuperLearner(Y=as.numeric(disc.U==bin), X=W, newX=newW, family='binomial', SL.library = SL.library, method='method.NNLS',
                                                 #potentially replace with
                                                 control = list(saveFitLibrary=saveFitLibrary), cvControl = list(V=n.folds, validRows=validRows)), silent=TRUE))
    }
    if(class(bin.fit) == "try-error") {
      capture.output(bin.fit <- try(SuperLearner(Y=as.numeric(disc.U==bin), X=W, newX=newW, family='binomial', SL.library = SL.library, method='method.NNLS2', control = list(saveFitLibrary=saveFitLibrary), cvControl = list(V=n.folds, validRows=validRows)), silent=TRUE))
    }
    if(class(bin.fit) != "try-error") {
      bin.fits[[paste0("bin", bin, ".SL")]] <- bin.fit
      bin.probs[,bin] <- bin.fit$SL.predict
    } else {
      bin.mean <- mean(as.numeric(disc.U==bin))
      if(class(SL.library) == "character") n.algs <- length(SL.library)
      else n.algs <- sum(unlist(lapply(SL.library, function(sl) length(sl) - 1)))
      bin.fits[[paste0("bin", bin, ".SL")]] <- list(Z = matrix(bin.mean, nrow=n,ncol=n.algs), library.predict = matrix(bin.mean, nrow=m,ncol=n.algs))
      bin.probs[,bin] <- bin.mean
    }
  }



  #SL.bin.probs <- .make.doubly.stochastic(bin.probs, row.sums = rep(1, m), col.sums = bin.fracs * m)
  SL.bin.probs <- bin.probs / rowSums(bin.probs)

  SL.densities <- SL.bin.probs[cbind(1:m, disc.U.new)] / bin.sizes[disc.U.new]

  n.alg <- ncol(bin.fits[["bin1.SL"]]$Z)
  cv.library.densities <-  matrix(NA, nrow=n, ncol=n.alg)
  library.densities <- matrix(NA, nrow=m, ncol=n.alg)
  for (j in 1:n.alg) {
    cv.bin.probs <- library.bin.probs <- matrix(NA, nrow=n, ncol=b)
    for (bin in 1:b) {
      cv.bin.probs[, bin] <- bin.fits[[paste0("bin", bin, ".SL")]]$Z[,j]
      library.bin.probs[, bin] <- bin.fits[[paste0("bin", bin, ".SL")]]$library.predict[,j]
    }
    if(any(is.na(cv.bin.probs)) | any(colSums(cv.bin.probs) == 0) | any(rowSums(cv.bin.probs) == 0)) {
      cv.library.densities[,j] <- rep(NA, n)
    } else {
      #cv.bin.probs <- .make.doubly.stochastic(cv.bin.probs, row.sums = rep(1, n), col.sums = bin.fracs * n)
      cv.bin.probs <- cv.bin.probs / rowSums(cv.bin.probs)
      cv.library.densities[,j] <- cv.bin.probs[cbind(1:n, disc.U)] / bin.sizes[disc.U]
    }

    if(any(is.na(library.bin.probs)) | any(colSums(library.bin.probs) == 0) | any(rowSums(library.bin.probs) == 0)) {
      library.densities[,j] <- rep(NA, m)
    } else {
      #library.bin.probs <- .make.doubly.stochastic(library.bin.probs, row.sums = rep(1, n), col.sums = bin.fracs * n)
      library.bin.probs <- library.bin.probs / rowSums(library.bin.probs)
      library.densities[,j] <- library.bin.probs[cbind(1:m, disc.U.new)] / bin.sizes[disc.U.new]
    }

  }

  alg.names <- paste0(bin.fits[["bin1.SL"]]$libraryNames, "_", b, "bins")

  ret <- list(bins = bins,  a.ecdf = a.ecdf, SL.bin.probs = SL.bin.probs, SL.densities = SL.densities, cv.library.densities = cv.library.densities, library.densities = library.densities, alg.names = alg.names)
  if(saveFitLibrary) ret$bin.fits <- bin.fits
  return(ret)

}

.make.doubly.stochastic <- function(mat, row.sums, col.sums, tol = .001) {
  ret <- mat
  while(sum(abs(rowSums(ret) - row.sums)) > tol | sum(abs(colSums(ret) - col.sums)) > tol) {
    ret <- ret / (rowSums(ret) / row.sums)
    ret <- t( t(ret) / (colSums(ret) / col.sums))
  }
  return(ret)
}

#' Prediction method for cmdSuperLearner
#'
#' This function predicts standardized conditional density function values given a fitted model and new data,
#'
#' @param fit Fitted \code{\link{cmdSuperLearner}} object. Must have been run with \code{control$saveFitLibrary = TRUE}.
#' @param newA \code{m x 1} numeric vector of new exposure values.
#' @param newW \code{m x p} data.frame of new covariate values.
#' @param threshold Minimum coefficient value for which library algorithms should be included in the prediction.
#' @return \code{cmdSuperLearner} returns a named list with the following elements:
#' \item{fit.times}{The time points at which the counterfactual survival curves (and contrasts) were fit.}
#'
#' @examples
#' # Sample data
#' set.seed(220)
#' # n <- 1000
#' # W <- data.frame(W1 = runif(n))
#' # Z <- rbinom(n, size = 1, prob = 1/(1 + exp(2-W$W1)))
#' # A <- (1-Z) * rnorm(n, mean = W$W1, sd = abs(1 + W$W1))
#' # fit <- cmdSuperLearner(A, W, control=list(SL.library = c("SL.mean", "SL.glm", "SL.gam", "SL.earth"), verbose=TRUE, n.bins = c(2:10)))
#' # Get predicted standardized "density" (really mass) at 0
#' # pred <- predict.cmdSuperLearner(fit = fit, newA = rep(0, n), newW = W)
#' # true.g <- 1 /(integrate(function(x) 1/(1 + exp(2-x)), 0, 1)$value * (1 + exp(2-W$W1)))
#' # plot(true.g, pred)
#' # abline(0,1, col='red')


predict.cmdSuperLearner <- function(fit, newA, newW, threshold = .001) {
  library(SuperLearner)
  newW <- as.data.frame(newW)
  new.U <- fit$a.ecdf(newA)
  trunc.coef <- fit$coef
  trunc.coef[trunc.coef < threshold] <- 0
  trunc.coef <- trunc.coef / sum(trunc.coef)
  nonzero <- which(trunc.coef > 0)
  lib.name.splits <- strsplit(fit$library.names, "_")
  lib.name.nbins <- unlist(lapply(lib.name.splits, function(l) as.numeric(strsplit(l[3], "bins")[[1]])))
  lib.name.alg <- unlist(lapply(lib.name.splits, function(l) paste0(l[1:2], collapse="_")))
  bins.to.fit <- unique(lib.name.nbins[nonzero])
  pred.densities <- matrix(NA, nrow=length(newA), ncol=length(fit$library.names))
  for(bin in bins.to.fit) {
    ind <- which(fit$control$n.bins == bin)
    new.bins <- .find.bin(new.U, bins = fit$fits[[ind]]$bins)
    bin.sizes <- unlist(lapply(fit$fits[[ind]]$bins, interval_measure))
    pred.probs <- matrix(NA, nrow = length(new.U), ncol = length(unique(lib.name.alg)))
    for(k in 1:length(fit$fits[[ind]]$bin.fits)) {
      if(any(new.bins == k)) {
        pred.probs[new.bins == k,] <- predict.SuperLearner(fit$fits[[ind]]$bin.fits[[k]], newdata = newW[new.bins == k,, drop=FALSE], onlySL = TRUE)$library.predict
      }
    }
    pred.probs <- pred.probs / rowSums(pred.probs)
    pred.densities[,which(lib.name.nbins == bin)] <- pred.probs / bin.sizes[new.bins]
  }
  c(pred.densities[,nonzero,drop=FALSE] %*% trunc.coef[nonzero])
}

#' Functions for converting from an IC and estimate of a function to the
#' corresponding basis function and IC.
#' @param est estimate of function
#' @param ic IC estimate
#'
#' @importFrom fda create.fourier.basis eval.basis
#' @export

convert_to_basis <- function(est, ic, gridvals, num_basis = 10) {
  uncor_est <- est - colMeans(ic)
  newfun <- stepfun(x = gridvals, y = c(0, uncor_est))
  fbasis_obj <- fda::create.fourier.basis(rangeval = c(0, 1),
                                          nbasis = num_basis)
  newgrid <- seq(0, 1, length.out = 1000)
  om_grid <- newfun(newgrid)
  mat <- fda::eval.basis(newgrid, fbasis_obj)
  ic_at_newgrid <- t(apply(ic, MARGIN = 1, FUN = function(x) {
    obs_fun <- stepfun(gridvals, y = c(0, x))
    return(obs_fun(newgrid))
  }))
  numr <- ic_at_newgrid %*% mat / nrow(mat)
  denom <- colMeans(mat ** 2)
  new_ic <- sweep(x = numr, MARGIN = 2, denom, "/")
  newest <- apply(mat, 2, function(x) mean(x * om_grid)) / denom
  return(list(est = as.numeric(newest + colMeans(new_ic)), ic = new_ic))
}
