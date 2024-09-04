###################################################################################################################################################
# Name: init_UGMM
#
## Description
#
## Arguments:
## X:          (n.obs x p) complete (numeric) data matrix
## G:          Number of components/clusters (numeric)
## m:          Number of variable groups (numeric)
## model.name: "FFFF" (string)
## initG:      Initialization for the unit partition (string)
## initm:      Initialization for the variable partition (string) 
##
## Values:
## label:      (n.obs x 1) vector of cluster labels
## model:      Trained model structure
## Sg:         List of G sample covariance matrices 
##################################################################################################################################################
init_UGMM <-
  function(X,
           G,
           m,
           model.name,
           initG,
           initm) {
    n.obs <- dim(X)[1]
    p <- dim(X)[2]
    Sigma <- vector(mode = "list", length = G)
    model <- list()
    ################################
    ##  POSTERIOR PROBABILITY     ##
    ################################
    if (initG == "kmeans") {
      w <- mclust::unmap(kmeans(X, G, iter.max = 100)$cluster)
    } else if (initG == "kmeansf") {
      w <- unname(ppclust::fcm(X, G, iter.max = 100, con.val = sqrt(.Machine$double.eps))$u)
    } else if (initG == "random"){
      w <- rand.member(n.obs, G)
    }
    w[which(w < .Machine$double.xmin, arr.ind = TRUE)] <- .Machine$double.xmin
    label <- apply(w, 1, which.max)
    ################################
    ##  VARIABLE PARTITION        ##
    ################################
    V <- vector(mode = "list", length = G)
    if (m == 1) {
      V[seq(1:G)] <- list(as.matrix(rep(1, p)))
    } else if (m == p) {
      V[seq(1:G)] <- list(diag(p))
    } else {
      if (initm == "ucms") {
        for (g in 1:G) {
          if (nrow(X[label == g, , drop = FALSE]) == 1) {
            V[[g]] <- rand.member(p, m)
          } else {
            S <- matrix(0, p, p)
            S <- cov(X[label == g, ]) 
            V[[g]] <- UCMSigmaV(S, m, rndst = 5)
          }
        }
      } else if  (initm == "random"){
        V[seq(1:G)] <- replicate(G, list(rand.member(p, m)))
      }
    }
    ################################
    ##    PRIOR PROBABILITIES     ##
    ################################
    model$pp <- prior(w)
    ################################################
    ##   COMPONENT MEAN VECTOR and SIGMA MATRIX   ##
    ################################################
    comp.param <- component_param(X, w)
    model$mu <- comp.param$mu
    Sigma <- comp.param$Sigma
    ################################################
    ##        COMPONENT ULTRAMETRIC SIGMA         ##
    ################################################
    for (g in 1:G) {
      Sv <- estimate_Sv(Sigma[[g]], V[[g]], model.name = model.name)
      Sw <- estimate_Sw(Sigma[[g]], Sv, V[[g]], model.name = model.name)
      Sb <- estimate_Sb(Sigma[[g]], V[[g]], model.name = model.name)
      if (m > 1) {
        Sw <- check_constraint_SwSb(Sw, Sb, model.name = model.name)$Sw
        if (m == p) {
          Sv <- Sw
        }
      }
      if (m < p) {
        Sv <- check_constraint_SvSw(Sv, Sw, model.name = model.name)$Sv
      }
      St <- pd_sigma(V[[g]] %*% (Sw + Sb) %*% t(V[[g]]) - diag(diag(V[[g]] %*% Sw %*% t(V[[g]]))) + diag(diag(V[[g]] %*% Sv %*% t(V[[g]]))))
      Sv <- Sv + St$a * diag(m)
      St <- St$S
      if (any(apply(V[[g]], 2, sum) == 1)) {
        Sw[which(apply(V[[g]], 2, sum) == 1), which(apply(V[[g]], 2, sum) == 1)] <- Sv[which(apply(V[[g]], 2, sum) == 1), which(apply(V[[g]], 2, sum) == 1)]
      }
      model$Sigma[[g]] <- St
      model$V[[g]] <- V[[g]]
      model$Sv[[g]] <- Sv
      model$Sw[[g]] <- Sw
      model$Sb[[g]] <- Sb
    }
    model$post <- w
    return(list(
      label = label,
      model = model,
      Sg = Sigma
    )
    )
  }
