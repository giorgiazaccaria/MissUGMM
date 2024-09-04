## Internal functions MissUGMM #################################################

tr <-
  function(x)
    sum(diag(as.matrix(x)))

################################################################################

norm <-
  function(x,
           type = c("none", "standard", "center", "range", "SVD")) {
    type <- match.arg(type,
                      choices = eval(formals(norm)$type),
                      several.ok = FALSE)
    x <- as.matrix(x)
    switch(
      type,
      "none" = x,
      "standard" = scale(x, center = apply(x, 2, mean, na.rm = TRUE), scale = apply(x, 2, sd, na.rm = TRUE)),
      "center"   = scale(x, center = apply(x, 2, mean, na.rm = TRUE), scale = FALSE),
      "range"    = apply(x, 2, function(x)
        (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))),
      "SVD"      = {
        x <- scale(x, center = apply(x, 2, mean, na.rm = TRUE), scale = apply(x, 2, sd, na.rm = TRUE))
        p <- ncol(x)
        SVD <- svd(x, nu = 0)
        x %*% SVD$v %*% diag(1 / sqrt(SVD$d), p, p)
      }
    )
  }

################################################################################

rand.member <-
  function(n.obs,
           G) {
    if (G > n.obs)
      stop('The number of groups is larger than the number of observations')
    if (G == n.obs) {
      U <- diag(n.obs)
    } else if (G == 1) {
      U <- c(rep(1, n.obs))
    }
    else {
      U <- matrix(0, n.obs, G)
      U[1:G,] = diag(G)
      U[(G + 1):n.obs, 1] <- 1
      for (i in (G + 1):n.obs) {
        U[i,] <- U[i, sample(G)]
      }
      U <- U[sample(n.obs),]
    }
    return(as.matrix(U))
  }

################################################################################

UCMSigmaV <-
  function(S,
           m,
           rndst,
           maxiter = 100,
           eps = sqrt(.Machine$double.eps),
           model.name = "FFFF") {
    p <- dim(S)[1]
    Im <- diag(m)
    
    if (m == p) {
      Vopt <- diag(p)
    } else if (m == 1) {
      Vopt <- as.matrix(rep(1, p))
    } else {
      for (loop in 1:rndst) {
        it <- 0
        V <- rand.member(p, m)
        Sv <- estimate_Sv(S, V, model.name)
        Sw <- estimate_Sw(S, Sv, V, model.name)
        Sb <- estimate_Sb(S, V, model.name)
        
        if (m > 1) {
          Sw <- check_constraint_SwSb(Sw, Sb, model.name)$Sw
          if (m == p) {
            Sv <- Sw
          }
        }
        if (m < p) {
          Sv <- check_constraint_SvSw(Sv, Sw, model.name)$Sv
        }
        
        St <- pd_sigma(V %*% (Sw + Sb) %*% t(V) - diag(diag(V %*% Sw %*% t(V))) + diag(diag(V %*% Sv %*% t(V))))
        Sv <- Sv + St$a * diag(m)
        St <- St$S
        
        fo <- tr(MASS::ginv(St,  tol = .Machine$double.xmin) %*% S) + log(det(St))
        fmin <- fo
        
        if (loop == 1) {
          Vopt <- V
          fopt <- fmin
        }
        
        fdif <- fmin
        while (fdif > eps && it <= maxiter) {
          it <- it + 1
          fo <- fmin
          for (i in 1:p) {
            posmin <- which(V[i, ] == 1)
            V[i,] <- Im[posmin, ]
            for (j in 1:m) {
              V[i,] <- Im[j, ]
              if (sum(V[, posmin]) > 0) {
                Sv <- estimate_Sv(S, V, model.name)
                Sw <- estimate_Sw(S, Sv, V, model.name)
                Sb <- estimate_Sb(S, V, model.name)
                
                if (m > 1) {
                  Sw <- check_constraint_SwSb(Sw, Sb, model.name)$Sw
                  if (m == p) {
                    Sv <- Sw
                  }
                }
                if (m < p) {
                  Sv <- check_constraint_SvSw(Sv, Sw, model.name)$Sv
                }
                
                St <- pd_sigma(V %*% (Sw + Sb) %*% t(V) - diag(diag(V %*% Sw %*% t(V))) + diag(diag(V %*% Sv %*% t(V))))
                Sv <- Sv + St$a * diag(m)
                St <- St$S
                
                ff <- tr(MASS::ginv(St,  tol = .Machine$double.xmin) %*% S) + log(det(St))
                
                if (ff < fmin) {
                  fmin <- ff
                  posmin <- j
                }
              }
            }
            V[i,] <- Im[posmin, ]
          }
          
          fdif <- fo - fmin
          if (fmin < fopt) {
            Vopt <- V
            fopt <- fmin
          }
        }
      }
    }
    return(Vopt)
  }

################################################################################

ultrcov <-
  function(Sb,
           V) {
    m <- nrow(Sb)
    A <- matrix(0, m - 1, 1)
    B <- matrix(0, m - 1, 1)
    levfus <- matrix(0, m - 1, 1)
    PP <- matrix(0, m, 1)
    P <- matrix(0, m, 1)
    class <- matrix(0, m, 1)
    SU <- matrix(0, m, m)
    
    for (q in 1:m) {
      PP[q, 1] <- q
      class[q, 1] <- sum(V[, q])
    }
    
    for (ustep in 1:(m - 1)) {
      rmax <- (-.Machine$double.xmax)
      
      for (i in 1:(m - 1)) {
        if (PP[i, 1] == i) {
          for (j in (i + 1):m) {
            if (PP[j, 1] == j) {
              if (Sb[i, j] >= rmax) {
                ic <- i
                jc <- j
                rmax <- Sb[i, j]
              }
            }
          }
        }
      }
      
      for (j in 1:m) {
        if (PP[j, 1] == jc) {
          PP[j, 1] <- ic
        }
      }
      
      A[ustep, 1] <- ic
      B[ustep, 1] <- jc
      levfus[ustep, 1] <- rmax
      
      for (i in 1:m) {
        if (i != ic && PP[i, 1] == i) {
          rs <-
            (class[ic, 1] * Sb[ic, i] + class[jc, 1] * Sb[jc, i]) / (class[ic, 1] + class[jc, 1])
          Sb[ic, i] <- rs
          Sb[i, ic] <- rs
        }
      }
      
      class[ic, 1] <- class[ic, 1] + class[jc, 1]
    }
    
    for (i in 1:m) {
      P[i, 1] <- i
    }
    
    for (k in 1:(m - 1)) {
      for (i in A[k, 1]:m) {
        if (P[i, 1] == A[k, 1]) {
          for (j in B[k, 1]:m) {
            if (P[j, 1] == B[k, 1]) {
              SU[i, j] <- levfus[k, 1]
              SU[j, i] <- levfus[k, 1]
            }
          }
        }
      }
      
      for (i in A[k, 1]:m) {
        if (P[i, 1] == B[k, 1]) {
          P[i, 1] <- A[k, 1]
        }
      }
    }
    return(SU)
  }

################################################################################

prior <-
  function(w) {
    pp <- colSums(w) / dim(w)[1]
    return(pp)
  }

################################################################################

post_missing <-
  function(X,
           difpat,
           pat,
           pp,
           mu,
           Sigma) {
    G <- dim(mu)[1]
    w <- matrix(as.double(NA), dim(X)[1], G)
    for (t in 1:nrow(difpat)){
      var.obs <- !difpat[t, , drop = FALSE]
      for (g in 1:G) {
        mu_o <- mu[g, var.obs, drop = FALSE]
        sigma_oo <- Sigma[[g]][var.obs, var.obs, drop = FALSE]
        w[pat == t, g] <-  log(pp[g]) + mclust::dmvnorm(X[pat == t, var.obs, drop = FALSE], mu_o, sigma_oo, log = TRUE)
      }
    }
    wnorm <- apply(w, 1, max)
    w <- exp(sweep(w, 1, wnorm, "-"))
    w <- w / rowSums(w)
    w[which(w < sqrt(.Machine$double.eps), arr.ind = TRUE)] <- sqrt(.Machine$double.eps)
    return(w)
  }

################################################################################

component_param <-
  function(X,
           w) {
    n.obs <- dim(X)[1]
    p <- dim(X)[2]
    G <- dim(w)[2]
    Sigma.new <- list()
    
    mu.new <- sweep(t(w) %*% X, 1, 1 / colSums(w), "*")
    
    r <- sqrt(w)
    for (g in 1:G) {
      Xo <- sweep(X, 2, mu.new[g,], "-", check.margin = FALSE)
      Xo <- sweep(Xo, 1, r[, g], "*", check.margin = FALSE)
      Sigma.new[[g]] <-  (t(Xo) %*% Xo) / colSums(w)[g]
    }
    return(list(mu = mu.new,
                Sigma = Sigma.new))
  }

################################################################################

impute_mis <-
  function(X,
           difpat,
           pat,
           pp,
           mu,
           Sigma) {
    n.obs <- dim(X)[1]
    p <- dim(X)[2]
    G <- dim(mu)[1]
    Ximp <- array(rep(X, G), dim = c(n.obs, p, G))
    sigma.imp <- array(0, dim = c(p, p, n.obs, G))
    for (t in 1:nrow(difpat)) {
      var.mis <- difpat[t, , drop = FALSE] != 0
      for (g in 1:G) {
        if (any(var.mis)) {
          var.obs <- !var.mis
          mu_m <- mu[g, var.mis, drop = FALSE]
          mu_o <- mu[g, var.obs, drop = FALSE]
          sigma_mo <- Sigma[[g]][var.mis, var.obs, drop = FALSE]
          sigma_om <- Sigma[[g]][var.obs, var.mis, drop = FALSE]
          sigma_mm <- Sigma[[g]][var.mis, var.mis, drop = FALSE]
          sigma_oo_inv <- MASS::ginv(Sigma[[g]][var.obs, var.obs], tol = .Machine$double.xmin)
          Ximp[pat == t, var.mis, g] <- 
            sweep(t(sigma_mo %*% sigma_oo_inv %*% t(sweep(X[pat == t, var.obs, drop = FALSE], 2, mu_o, "-", check.margin = FALSE))), 2, mu_m, "+")
          sigma.mis <- sigma_mm - sigma_mo %*% sigma_oo_inv %*% sigma_om
          un.t <- which(pat==t)
          for (i in 1:length(un.t)) {
            sigma.imp[var.obs, var.obs, un.t[i], g] <- tcrossprod(t(X[un.t[i], var.obs, drop = F] - mu_o))
            sigma.imp[var.obs, var.mis, un.t[i], g] <- tcrossprod(t(X[un.t[i], var.obs, drop = F] - mu_o), t(t(as.matrix(Ximp[un.t[i], var.mis, g, drop = FALSE])) - mu_m))
            sigma.imp[var.mis, var.obs, un.t[i], g] <- t(sigma.imp[var.obs, var.mis, un.t[i], g])
            sigma.imp[var.mis, var.mis, un.t[i], g] <- tcrossprod(t(t(as.matrix(Ximp[un.t[i], var.mis, g, drop = FALSE])) - mu_m)) + sigma.mis
          }
        } else {
          un.t <- which(pat==t)
          for (i in 1:length(un.t)) {
            sigma.imp[, , un.t[i], g] <- tcrossprod(t(X[un.t[i], , drop = FALSE] - mu[g, , drop = FALSE]))
          }
        }
      }
    }
    return(list(
      Ximp = Ximp,
      sigma.imp = sigma.imp))
  }

################################################################################

component_param_missing <-
  function(Ximp,
           w, 
           sigma.imp) {
    n.obs <- dim(Ximp)[1]
    p <- dim(Ximp)[2]
    G <- dim(w)[2]
    mu.new <- matrix(0, G, p)
    Sigma.new <- list()
    
    for (g in 1:G){
      mu.new[g, ] <- sweep(t(w[, g]) %*% Ximp[, , g], 2, 1 / sum(w[, g]), "*")
      Xo <- apply(sweep(sigma.imp[, , , g], 3, w[, g], "*", check.margin = FALSE), 1:2, sum)
      Sigma.new[[g]] <-  Xo / sum(w[, g])
    }
    return(list(mu = mu.new,
                Sigma = Sigma.new))
  }

################################################################################

loglik_EI_G_pugmm <-
  function(X,
           pp,
           mu,
           Sigma,
           Sv,
           Sw,
           Sb,
           V,
           gaussian) {
    n.obs <- dim(X)[1]
    p <- dim(X)[2]
    G <- length(pp)
    llh <- rep(as.double(NA), G)
    lf <- matrix(as.double(NA), n.obs, G)
    for (g in 1:G) {
      nv <- colSums(V[[g]])
      if (any(nv == 1)) {
        gaussian = "mclust"
      }
      if (gaussian == "mclust") {
        lf[, g] <-
          mclust::dmvnorm(X, mu[g, , drop = FALSE], Sigma[[g]], log = TRUE) + log(pp[g])
      } else {
        lf[, g] <-
          canonical_norm(X, mu[g, , drop = FALSE], Sv[[g]], Sw[[g]], Sb[[g]], V[[g]]) + log(pp[g])
      }
    }
    loglik <- sum(log(rowSums(exp(lf))))
    return(loglik)
  }

################################################################################

loglik_F_pugmm <-
  function(X,
           g,
           pp,
           mu,
           Sigma.current,
           Sigma.other,
           Sv.current,
           Sv.other,
           Sw.current,
           Sw.other,
           Sb.current,
           Sb.other,
           V.current,
           V.other,
           gaussian) {
    n.obs <- dim(X)[1]
    G <- length(pp)
    llh <- rep(as.double(NA), G)
    lf <- matrix(as.double(NA), n.obs, G)
    nv0 <- colSums(V.current)
    for (k in 1:G) {
      nv <- colSums(V.other[[k]])
      if (any(nv0 == 1) || any(nv == 1)) {
        gaussian = "mclust"
      }
      if (gaussian == "mclust") {
        if (k == g) {
          lf[, k] <-
            as.matrix(mclust::dmvnorm(X, mu[k, , drop = FALSE], Sigma.current, log = TRUE)) + log(pp[k])
        } else {
          lf[, k] <-
            as.matrix(mclust::dmvnorm(X, mu[k, , drop = FALSE], Sigma.other[[k]], log = TRUE)) + log(pp[k])
        }
      } else {
        for (k in 1:G) {
          if (k == g) {
            lf[, k] <-
              canonical_norm(X, mu[k, , drop = FALSE], Sv.current, Sw.current, Sb.current, V.current) + log(pp[k])
          } else {
            lf[, k] <-
              canonical_norm(X, mu[k, , drop = FALSE], Sv.other[[k]], Sw.other[[k]], Sb.other[[k]], V.other[[k]]) + log(pp[k])
          }
        }
      }
    }
    loglik <- sum(log(rowSums(exp(lf))))
    return(loglik)
  }

################################################################################

loglik_MissUGMM <-
  function(X,
           difpat,
           pat,
           pp,
           mu,
           Sigma) {
    n.obs <- dim(X)[1]
    p <- dim(X)[2]
    G <- length(pp)
    lf <- matrix(as.double(NA), n.obs, G)
    for (t in 1:nrow(difpat)){
      var.obs <- !difpat[t, , drop = FALSE]
      for (g in 1:G) {
        mu_o <- mu[g, var.obs, drop = FALSE]
        sigma_oo <- Sigma[[g]][var.obs, var.obs, drop = FALSE]
        lf[pat == t, g] <-
          mclust::dmvnorm(X[pat == t, var.obs, drop = FALSE], mu_o, sigma_oo, log = TRUE) + log(pp[g])
      }
    }
    loglik <- sum(log(rowSums(exp(lf))))
    return(loglik)
  }

################################################################################

loglik_MissUGMM_int <-
  function(X,
           g,
           difpat,
           pat,
           pp,
           mu,
           Sigma.current,
           Sigma.other) {
    n.obs <- dim(X)[1]
    G <- length(pp)
    lf <- matrix(as.double(NA), n.obs, G)
    for (t in 1:nrow(difpat)){
      var.obs <- !difpat[t, , drop = FALSE]
      for (k in 1:G) {
        mu_o <- mu[k, var.obs, drop = FALSE]
        if (k == g){
          sigma_oo <- Sigma.current[var.obs, var.obs, drop = FALSE]
        } else {
          sigma_oo <- Sigma.other[[k]][var.obs, var.obs, drop = FALSE]
        }
        lf[pat == t, k] <-
          mclust::dmvnorm(X[pat == t, var.obs, drop = FALSE], mu_o, sigma_oo, log = TRUE) + log(pp[k])
      }
    }
    loglik <- sum(log(rowSums(exp(lf))))
    return(loglik)
  }

################################################################################

pd_sigma <-
  function(S) {
    p <- dim(S)[1]
    a <- as.numeric(0)
    while (min(eigen(S)$values) <= 0 ||
           det(S) <= .Machine$double.xmin) {
      a <- a + abs(min(eigen(S)$values)) + sqrt(.Machine$double.eps)
      S <- S + a * diag(p)
    }
    return(list(S = S + a * diag(p),
                a = a))
  }

################################################################################
aitken_obj <-
  function(vloglik) {
    if (length(vloglik) < 3) {
      return(Inf)
    } else {
      l.loglik <- length(vloglik)
      loglik.c <- vloglik[l.loglik]
      loglik <- vloglik[(l.loglik - 1)]
      loglik.p <- vloglik[(l.loglik - 2)]
      ait.c <- (loglik.c - loglik) / (loglik - loglik.p)
      loglik.inf <- loglik  + ((loglik.c - loglik)/(1 - ait.c))
      conv <- loglik.inf - loglik
      if (conv < 0) {
        conv <- 1
      }
      return(conv)
    }
  }

################################################################################

number_param <-
  function(p,
           G,
           m,
           model.name = c("EUUU", "EUUE", "EUEE", "EEEU", "EEEE", "EEEF", "EEFF", "EFFF", "FIII", "FIIF", "FIFF", "FFFI", "FFFF")) {
    model.name <- match.arg(model.name)
    switch(
      model.name,
      "EUUU" = {pm <- G + (G + 1) * p + 2},
      "EUUE" = {pm <- G + (G + 1) * p + m},
      "EUEE" = {pm <- G + (G + 1) * p + 2* m - 1},
      "EEEU" = {pm <- G + ((G + 1) * p) + 2 * m},
      "EEEE" = {pm <- G + (G + 1) * p + 3 * m - 2},
      "EEEF" = {pm <- (G + 1) * p + (G + 2) * m - 1},
      "EEFF" = {pm <- (G + 1) * p + (2 * G  + 1) * m - 1},
      "EFFF" = {pm <- (G + 1) * p + 3 * G * m - 1},
      "FIII" = {pm <- 2 * G * (p + 2) - 1},
      "FIIF" = {pm <- G * (2 * p + m + 2) - 1},
      "FIFF" = {pm <- G * (2 * p + 2* m + 1) - 1},
      "FFFI" = {pm <- 2 * G * (p + m + 1) - 1},
      "FFFF" = {pm <- G * (2 * p + 3 * m) - 1}
    )
    return(pm)
  }



