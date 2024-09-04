source('Initialization_UGMM.R')
source('InternalFunctions_MissUGMM.R')
source('CovarianceStructures.R')

###################################################################################################################################################
# Name: MissUGMM
#
## Description
#
## Arguments:
## X:                 (n.obs x p) data matrix with missing values
## G:                 Number of components/clusters (numeric)
## m:                 Number of variable groups (numeric)
## normalization:     Type of normalization for the data : NULL, "standard", "center", "range" "SVD" (default: NULL, string)
## maxiter:           Maximum number of iterations (default: 500, numeric)
## tol:               Tolerance value  for convergence (default: sqrt(.Machine$double.eps), numeric)
## stop:              "relative" = relative log-likelihood in two sequential iterations; "aitken" = Aitken acceleration-based stopping rule (default: "aitken", string)
## rndstart:          number of random starts (default: 20, numeric)
## initX:             "knn" = imputation of the missing values via k-nearest neighbors algorithm with k = 5;  "complete cases" = initialization only on complete observations; 
#                     "overall mean" = imputation of the missing values via the mean of the observed values per variable; "cluster-wise mean" = imputation of the missing values via the mean of the observed values per variable and cluster (randomly generated)  (default: "knn", string)
## initG:             "kmeans" = hard k-means with 1 random start; "kmeansf" = fuzzy k-means with 1 random start; "random" = random initialization (default: "kmeans", string)
## initm:             "ucms" = UCMSigma with 5 random start; "random" = random initialization  (default: "ucms", string)
#
## Values:
## label:             (n.obs x 1) vector of cluster membership
## pp:                (1 x G) prior probability vector
## mu:                (G x p) matrix of the component mean vectors (row by row)
## sigma:             List of G (p x p) extended ultrametric covariance matrices
## V:                 List of G (p x m) variable-group membership matrices
## Sv:                List of G (m x m) diagonal matrices of the group variances
## Sw:                List of G (m x m) matrix of the within-group covariances 
## Sb:                List of G (m x m) matrix of the between-group covariances 
## post:              (n.obs x G) matrix of the posterior probabilities
## pm:                Number of the model parameters
## pm.cov:            Number of the covariance parameters
## pm.free:           Number of free parameters (pm - constraints)
## count.constr.SwSb: Overall number of times the constraint between Sw and Sb has been turned on
## count.constr.SvSw: Overall number of times the constraint between Sv and Sw has been turned on
## bic:               Bayesian information criterion
## loglik:            Final value of the likelihood
## loop:              Loop corresponding to the best model
## iter:              Actual number of iterations needed to reach convergence
##################################################################################################################################################
MissUGMM <-
  function(X,
           G,
           m,
           normalization = NULL,
           maxiter = 500,
           tol = sqrt(.Machine$double.eps),
           stop = "aitken",
           rndstart = 20,
           initX = "knn",
           initG = "kmeans",
           initm = "ucms",
           showprogress = TRUE) {
    call <- match.call()
    model.name = "FFFF"
    if (!is.null(normalization)) {
      X <- norm(X, normalization)
    }
    else {
      X <- as.matrix(X)
    }
    if (!any(is.na(X))) {
      stop("The data set has no missing values", call. = TRUE)
    }
    n.obs <- dim(X)[1]
    p <- dim(X)[2]
    if (G < 1 || G > n.obs) {
      stop("G is not properly fixed: G must be chosen in [1, dim(X)[1]]", call. = TRUE)
    }
    if (m < 1 || m > p) {
      stop("m is not properly fixed: m must be chosen in [1, p]", call. = TRUE)
    }
    IQ <- diag(m)
    wpc <- G * (G + 1) / 2
    # PRE-PROCESSING FOR THE INITIALIZATION      
    pattern <- matrix(as.numeric(is.na(X)), ncol = p)
    dif.pattern <- pattern[!duplicated(pattern), , drop = FALSE]
    if (any(apply(dif.pattern, 1, sum) == p)) {
      stop("The data set has one or more observations entirely missing", call. = TRUE)
    }
    pat.unit <- matrix(0, n.obs, 1)
    for (t in 1:nrow(dif.pattern)) {
      pat.unit[apply(pattern, 1, identical, dif.pattern[t, ]), ] <- t
    }
    if (initX == "knn") {
      Ximp <- bnstruct::knn.impute(X, k = 5, cat.var = NULL)
    } else if (initX == "complete cases") {
      Ximp <- X[complete.cases(X), ]
      if (dim(Ximp)[1] < p + 1) {
        stop("The number of complete cases is too small", call. = TRUE)
      }
    } else if (initX == "overall mean" || initX == "cluster-wise mean") {
      Ximp <- X
      mis.var <- apply(pattern, 2, function(x) return(any(x == 1)))
      if (initX == "overall mean") {
        Ximp[, mis.var] <- apply(X[, mis.var], 2, function(x) {
          x[which(is.na(x))] <- mean(x, na.rm = TRUE) 
          return(x)})
      } else {
        init.part <- rand.member(n.obs, G)
        for (g in 1:G) {
          Ximp[init.part[, g] == 1, mis.var] <- apply(X[init.part[, g] == 1, mis.var], 2, function(x) {
            x[which(is.na(x))] <- mean(x, na.rm = TRUE) 
            return(x)})
        }
      }
    }
    # STARTING
    for (loop in 1:rndstart) {
      count <- 1
      conv <- 1
      wp <- 0
      # INITIALIZATION  
      init <- init_UGMM(Ximp, G, m, model.name = model.name, initG = initG, initm = initm)
      model <- init$model
      label <- init$label
      Sigma <- init$Sg
      loglik <-
        loglik_MissUGMM(
          X,
          dif.pattern,
          pat.unit,
          model$pp,
          model$mu,
          model$Sigma)
      if (showprogress) {
        pb <- txtProgressBar(min = 0, max = rndstart, style = 3, 
                             width = 75, char = "=")
      }
      # ITERATIONS    
      while (conv > tol && count <= maxiter) {
        # E-STEP    
        model$post <- post_missing(X, dif.pattern, pat.unit, model$pp, model$mu, model$Sigma)
        label <- apply(model$post, 1, which.max)
        fmax <- loglik
        if (sum(unique(label)) != wpc) {
          wp <- 1
          loglik <- -.Machine$double.xmax
          break
        }
        val.impute <- impute_mis(X, dif.pattern, pat.unit, model$pp, model$mu, model$Sigma)
        model$X.imputed <- apply(rep(model$post, p) * aperm(val.impute$Ximp, perm = c(1, 3, 2)), 3, rowSums)
        # M-STEP   
        model$pp <- prior(model$post)
        comp.param <- component_param_missing(val.impute$Ximp, model$post, val.impute$sigma.imp)
        model$mu <- comp.param$mu
        Sigma <- comp.param$Sigma
        constr.SwSb <- 0 
        constr.SvSw <- 0 
        for (g in 1:G) {
          Stit <- model$Sigma[[g]]
          Svit <- model$Sv[[g]]
          Swit <- model$Sw[[g]]
          Sbit <- model$Sb[[g]]
          constr.SwSbit <- 0
          constr.SvSwit <- 0
          V0 <- model$V[[g]]
          for (i in 1:p) {
            posmax <- which(V0[i, ] == 1)
            V0[i,] <- IQ[posmax, ]
            for (h in 1:m) {
              count.SwSb <- 0
              count.SvSw <- 0
              V0[i, ] <- IQ[h, ]
              if (sum(V0[, posmax]) > 0) {
                Sv <- estimate_Sv(Sigma[[g]], V0, model.name)
                Sw <- estimate_Sw(Sigma[[g]], Sv, V0, model.name)
                Sb <- estimate_Sb(Sigma[[g]], V0, model.name)
                if (m > 1) {
                  Sw.constr <- check_constraint_SwSb(Sw, Sb, model.name)
                  Sw <- Sw.constr$Sw
                  count.SwSb <- Sw.constr$constr.SwSb
                  if (m == p) {
                    Sv <- Sw
                  }
                }
                if (m < p) {
                  Sv.constr <- check_constraint_SvSw(Sv, Sw, model.name)
                  Sv <- Sv.constr$Sv
                  count.SvSw <- Sv.constr$constr.SvSw
                }
                St <- pd_sigma(V0 %*% (Sw + Sb) %*% t(V0) - diag(diag(V0 %*% Sw %*% t(V0))) + diag(diag(V0 %*% Sv %*% t(V0))))
                Sv <- Sv + St$a * diag(m)
                St <- St$S
                ff <-
                  loglik_MissUGMM_int(
                    X,
                    g,
                    dif.pattern,
                    pat.unit,
                    model$pp,
                    model$mu,
                    St,
                    model$Sigma)
                if (ff >= fmax) {
                  fmax <- ff
                  posmax <- h
                  Stit <- St
                  Svit <- Sv
                  Swit <- Sw
                  Sbit <- Sb
                  constr.SwSbit <- count.SwSb
                  constr.SvSwit <- count.SvSw
                }
              }
            }
            V0[i,] <- IQ[posmax, ]
          }
          if (any(apply(V0, 2, sum) == 1)) {
            Swit[which(apply(V0, 2, sum) == 1), which(apply(V0, 2, sum) == 1)] <- Svit[which(apply(V0, 2, sum) == 1), which(apply(V0, 2, sum) == 1)]
          }
          model$Sigma[[g]] <- Stit
          model$V[[g]] <- V0
          model$Sv[[g]] <- Svit
          model$Sw[[g]] <- Swit
          model$Sb[[g]] <- Sbit
          constr.SwSb <- constr.SwSb + constr.SwSbit
          constr.SvSw <- constr.SvSw + constr.SvSwit
        }
        # OBJECTIVE FUNCTION 
        loglik.c <-
          loglik_MissUGMM(
            X,
            dif.pattern,
            pat.unit,
            model$pp,
            model$mu,
            model$Sigma)
        if (stop == "aitken") {
          if (count < 4) {
            loglik.p <- loglik
          } else {
            if (loglik.c > loglik) {
              ait.c <- (loglik.c - loglik) / (loglik - loglik.p)
            } else {
              ait.c <- 0
            }
            loglik.inf <- loglik  + (loglik.c - loglik)/(1 - ait.c)
            conv <- loglik.inf - loglik
            loglik.p <- loglik
          }
        } else if (stop == "relative") {
          conv <- (loglik.c - loglik) / abs(loglik)
        }
        count <- count + 1
        loglik <- loglik.c
      }
      if (loop == 1) {
        if (wp == 1) {
          iter.best <- count
          constr.SwSb <- 0
          constr.SvSw <- 0
        }
        label.best <- label
        model.best <- model
        loglik.best <- loglik
        constr.SwSb.best <- constr.SwSb
        constr.SvSw.best <- constr.SvSw
        loop.best <- loop
        iter.best <- count - 1
      }
      if (loglik > loglik.best) {
        label.best <- label
        model.best <- model
        loglik.best <- loglik
        constr.SwSb.best <- constr.SwSb
        constr.SvSw.best <- constr.SvSw
        loop.best <- loop
        iter.best <- count - 1
      }
      if (G == 1) {
        break
      }
      if (showprogress) {
        setTxtProgressBar(pb, loop)
        cat("Loop", loop, "/", rndstart)
      }
    }
    pm.best <- number_param(p, G, m, model.name)
    free.param <- pm.best - (G * m + constr.SwSb.best + constr.SvSw.best)
    bic.best <- 2 * loglik.best - (free.param * log(n.obs))
    if (loglik.best == -.Machine$double.xmax) {
      print("The solution has a number of cluster < G and the objective is not computed.")
    }
    return(
      list(
        call = call,
        X.imputed = model.best$X.imputed,
        G = G,
        m = m,
        label = label.best,
        pp = model.best$pp,
        mu = model.best$mu,
        sigma = model.best$Sigma,
        V = model.best$V,
        Sv = model.best$Sv,
        Sw = model.best$Sw,
        Sb = model.best$Sb,
        post = model.best$post,
        pm = pm.best,
        pm.cov = pm.best - (G - 1 + G * p),
        pm.free = free.param,
        count.constr.SwSb = constr.SwSb.best,
        count.constr.SvSw = constr.SvSw.best,
        bic = bic.best,
        loglik = loglik.best,
        loop = loop.best,
        iter = iter.best
      )
    )
  }
