##=== extending MBN using INLA for hirarchecal (multilevel) data (correlated random effects $m \times m$ random effects ====#
library(bnlearn)
library(igraph)
library(INLA)
library(gridExtra)
library(ggplot2)
library(lme4)
library(dplyr)
library(mvtnorm)
library(psych)
library(caret)
library(graphpcor)

nsims = 100
##==== Custom score function MLE ====##
fitME = function(node, parents, group, data) {
  
  rhs = paste(c("1", parents), collapse = " + ")
  formula = paste(node, "~", rhs, "+", paste0("(", rhs, " | ", group, ")"))
  
  model = try(lmer(formula, data = data, control = lmerControl(calc.derivs = FALSE)))
  
  if (is(model, "try-error"))
    return(NULL)
  
  return(model)
  
}#FITME

scoreME =  function(node, parents, data, args = list()) {
  
  if (node == args$group)
    return(0)
  
  # mixed effects models with performance tuning to speed up learning.
  require(lme4)
  environment(nloptwrap)$defaultControl$maxeval = 5 * 1e3
  environment(nloptwrap)$defaultControl$xtol_abs = 1e-6
  environment(nloptwrap)$defaultControl$ftol_abs = 1e-6
  
  local.distribution =
    fitME(node = node, parents = setdiff(parents, args$group),
          group = args$group, data = data)
  
  if (is.null(local.distribution))
    return(-Inf)
  
  return(- BIC(local.distribution))
  
}#SCOREME

##=== Custom score function INLA_ default prior ===#

fitME1 = function(node, parents, group, data) {
  nid <- length(unique(data[[group]]))  
  data$numid <- as.numeric(as.factor(data[[group]]))  
  
  if (length(parents) == 0) {
    formula = as.formula(paste(node, "~ 1 + f(numid, model = 'iidkd', n =", nid, ", order = 2)"))
  } else {
    
    rhs = paste("1", paste(parents, collapse = " + "), sep = " + ")  
    m = length(parents) + 1  
    
    slope_names <- paste0("slopeid", seq_along(parents)) 
    data[[slope_names[1]]] <- data$numid + nid  
    
    if (length(parents) > 1) {
      for (i in 2:length(parents)) {
        data[[slope_names[i]]] <- data[[slope_names[i - 1]]] + nid  
      }
    }
    
    slope_terms <- paste0("f(", slope_names, ", ", parents, ", copy = 'numid')", collapse = " + ")
    
    formula = as.formula(paste(node, "~", rhs, "+ f(numid, model = 'iidkd', n =", m * nid, 
                               ", order =", m, ") +", slope_terms))
  }
  model = tryCatch(
    inla(formula, family = "gaussian", data = data,
         control.inla = list(int.strategy = "eb")),
    error = function(e) {
      message("Error fitting model for node ", node, ": ", e$message)
      return(NULL)
    }
  )
  
  return(model)
}

scoreME1 =  function(node, parents, data, args = list()) {
  
  if (node == args$group)
    return(0)
  
  # mixed effects models with performance tuning to speed up learning.
  require(INLA)
  
  local.distribution =
    fitME1(node = node, parents = setdiff(parents, args$group),
           group = args$group, data = data)
  
  if (is.null(local.distribution))
    return(-Inf)
  
  return(local.distribution$mlik[1])
  
  
}

#== score function INLA for default LKJ prior ===# 
fitME2 <- function(node, parents, group, data) {
  data[[group]] <- as.numeric(data[[group]])
  data$rintercept <- rep(1, each = ncluster * csize)
  
  eta_values <- c(10, 20) 
  
  if (length(parents) == 0) {
    m <- 2
  } else {
    m <- length(parents) + 1
    #if (length(parents) > 2) eta_values <- c(10, 20)
  }
  
  cinla <- list(int.strategy = 'eb')
  cfam <- list(hyper = list())
  
  # Build the formula based on number of parents
  if (length(parents) == 0) {
    formula_text <- paste0(node, " ~ 1 + f(rintercept, model = cmodel, replicate = ", group, ")")
  } else {
    rhs <- paste("1", paste(parents, collapse = " + "), sep = " + ")
    slope_names <- paste0("rSlope", seq_along(parents))
    
    data[[slope_names[1]]] <- rep(2, each = ncluster * csize)
    if (length(parents) > 1) {
      for (i in 2:length(parents)) {
        data[[slope_names[i]]] <- rep(i + 1, each = ncluster * csize)
      }
    }
    slope_terms <- paste0("f(", slope_names, ", ", parents, ", copy = 'rintercept', replicate = ", group, ")", collapse = " + ")
    
    formula_text <- paste0(node, " ~ ", rhs,
                           " + f(rintercept, model = cmodel, replicate = ", group, ")",
                           " + ", slope_terms)
  }
  
  for (eta in eta_values) {
    #message("Trying eta = ", eta, " for node ", node)
    
    cmodel <- cgeneric(
      model = "LKJ",
      n = m,
      eta = eta
    )
    
    formula <- as.formula(formula_text)
    model <- try(
      {
        inla(formula, family = "gaussian", data = data,
             control.inla = cinla,
             control.family = cfam)
      },
      silent = TRUE
    )
    
    if (!inherits(model, "try-error")) {
      return(model)
    } else {
      #message("INLA crashed for node ", node, " with eta = ", eta, ". Trying next eta...")
    }
  }
  
  warning("All model attempts failed for node ", node)
  return(NULL)
}


scoreME2 =  function(node, parents, data, args = list()) {
  
  if (node == args$group)
    return(0)
  require(INLA)
  
  local.distribution =
    fitME2(node = node, parents = setdiff(parents, args$group),
           group = args$group, data = data)
  
  if (is.null(local.distribution))
    return(-Inf)
  
  return(local.distribution$mlik[1])
}


#====Simulation study for m by m correlated RE with different Wishart priors=====#

#== Random DAG used to simulate the data ===#
TMBN = model2network("[G][Y1|G][Y2|G][Y8|G][Y3|Y1:G][Y4|Y1:G][Y5|Y2:G][Y6|Y1:Y5:G][Y9|Y1:Y4:Y8:G][Y7|Y5:Y6:G][Y10|Y4:Y6:G]")

shd_lme = list() 
shd_inla = list() 
shd_inla1 = list()  
 
for(k in 1:nsims) {
  ncluster = 30
  csize = 5
  N = ncluster*csize
  data = data.frame("ID" = rep(1:N,
                               each = 1),
                    "G" = rep(1:ncluster,
                              each = csize))
  #== Y1|G
  
  
  sigmau01 = 0.35 
  u0 = rnorm(n = ncluster,
              mean = 0,
              sd = sigmau01)
  
  data$u0 = rep(u0,
                 each = 5, times=1)
  
  sigmae = sqrt(0.5)
  data$e = rnorm(N, 0, sigmae)
  
  data$b0 = 0.72
  data$Y1 = data$b0 +  data$u0 + data$e
  
  data <- data[, !(names(data) %in% c("u0","u1", "e", "b0"))]
  
  #== Y2|Y1:G
  
  m= 2
  rho <- -0.4
  Sigma <- matrix(NA, m, m)
  diag(Sigma) <- c(0.4,0.35)
  for(i in 1:m) {
    for (j in 1:m) {
      if (i != j) {
        Sigma[i, j] <- rho^abs(i-j) * sqrt(Sigma[i, i] * Sigma[j, j])
      }
    }
  }
  random_effects <- rmvnorm(n = ncluster,sigma = Sigma)
  
  data$u0 = rep(random_effects[, 1],
                each = 5, times=1)
  
  data$u1 = rep(random_effects[, 2],
                each = 5, times=1)
  
  sigmae = sqrt(0.5) # 0.6
  data$e = rnorm(N, 0, sigmae)
  
  data$b0 = 1.3
  data$b1 = 0.3
  
  data$Y2 = data$b0 + data$b1*data$Y1 + data$u1*data$Y1 + data$u0 + data$e
  
  data <- data[, !(names(data) %in% c("u0","u1", "e", "b0", "b1"))]
  
  #== Y3|G
  sigmau0 = 0.2 # 0.22
  u0 = rnorm(n = ncluster,
              mean = 0,
              sd = sigmau0)
  
  data$u0 = rep(u0,
                 each = 5, times=1)
  
  sigmae = sqrt(0.52)
  data$e = rnorm(N, 0, sigmae)
  
  data$b0 = -1.03
  data$Y3 = data$b0 +  data$u0 + data$e
  
  data <- data[, !(names(data) %in% c("u0","u1", "e", "b0"))]
  
  #===Y4|Y2:G
  m= 2
  rho <- 0.54
  Sigma <- matrix(NA, m, m)
  diag(Sigma) <- c(0.48, 0.4)
  for(i in 1:m) {
    for (j in 1:m) {
      if (i != j) {
        Sigma[i, j] <- rho^abs(i-j) * sqrt(Sigma[i, i] * Sigma[j, j])
      }
    }
  }
  random_effects <- rmvnorm(n = ncluster,  sigma = Sigma)
  
  
  data$u0 = rep(random_effects[, 1],
                each = 5, times=1)
  
  data$u1 = rep(random_effects[, 2],
                each = 5, times=1)
  
  sigmae = sqrt(0.54) # 0.35
  data$e = rnorm(N, 0, sigmae)
  
  data$b0 = -0.32
  data$b1 = 0.56
  
  data$Y4 = data$b0 + data$b1*data$Y2 + data$u1*data$Y2 + data$u0 + data$e
  
  data <- data[, !(names(data) %in% c("u0","u1", "e", "b0", "b1"))]
  
  #== Y5|Y3:G 
  m= 2
  rho <- 0.4
  Sigma <- matrix(NA, m, m)
  diag(Sigma) <- c(0.52, 0.4)
  for(i in 1:m) {
    for (j in 1:m) {
      if (i != j) {
        Sigma[i, j] <- rho^abs(i-j) * sqrt(Sigma[i, i] * Sigma[j, j])
      }
    }
  }
  random_effects <- rmvnorm(n = ncluster,  sigma = Sigma)
  
  
  
  data$u0 = rep(random_effects[, 1],
                each = 5, times=1)
  
  data$u1 = rep(random_effects[, 2],
                each = 5, times=1)
  
  sigmae = sqrt(0.450) #0.95
  data$e = rnorm(N, 0, sigmae)
  
  data$b0 = 2.63
  data$b1 = 1.19
  
  data$Y5 = data$b0 + data$b1*data$Y3 + data$u1*data$Y3 + data$u0 + data$e
  
  data <- data[, !(names(data) %in% c("u0","u1", "e", "b0", "b1i", "b1"))]
  
  #== Y6|Y1:Y3:G
  m= 3
  rho <- 0.3
  Sigma <-  matrix(NA, m, m)
  diag(Sigma) <- c(0.54, 0.43, 0.34)
  for(i in 1:m) {
    for (j in 1:m) {
      if (i != j) {
        Sigma[i, j] <- rho^abs(i-j) * sqrt(Sigma[i, i] * Sigma[j, j])
      }
    }
  }
  random_effects <- rmvnorm(n = ncluster,  sigma = Sigma)
  
  data$u0 = rep(random_effects[, 1],
                each = 5, times=1)
  
  data$u1 = rep(random_effects[, 2],
                each = 5, times=1)
  data$u2 = rep(random_effects[, 3],
                each = 5, times=1)
  
  sigmae = sqrt(0.4) #0.4
  data$e = rnorm(N, 0, sigmae)
  
  data$b0 = 1.03
  data$b1 = -0.27
  data$b2 = 0.77
  
  data$Y6 = data$b0 + data$b1*data$Y1 +data$b2*data$Y3 + data$u1*data$Y1 + 
    data$u2*data$Y3 + data$u0 + data$e
  
  data <- data[, !(names(data) %in% c("u0","u1","u2", "e", "b0", "b1", "b2"))]
  
  #=== Y7|Y4:G
  
  m= 2
  rho <- 0.4
  Sigma <- matrix(NA, m, m)
  diag(Sigma) <- c(0.54, 0.41)
  for(i in 1:m) {
    for (j in 1:m) {
      if (i != j) {
        Sigma[i, j] <- rho^abs(i-j) * sqrt(Sigma[i, i] * Sigma[j, j])
      }
    }
  }
  random_effects <- rmvnorm(n = ncluster,  sigma = Sigma)
  
  data$u0 = rep(random_effects[, 1],
                each = 5, times=1)
  
  data$u1 = rep(random_effects[, 2],
                each = 5, times=1)
  
  sigmae = sqrt(0.45) # 0.6
  data$e = rnorm(N, 0, sigmae)
  
  data$b0 = 2.2
  data$b1 = 0.98
  
  data$Y7 = data$b0 + data$b1*data$Y4 + data$u1*data$Y4 + data$u0 + data$e
  
  data <- data[, !(names(data) %in% c("u0","u1", "e", "b0", "b1"))]
  
  #== Y8|Y5:G
  
  m= 2
  rho <- 0.52
  Sigma <- matrix(NA, m, m)
  diag(Sigma) <- c(0.45, 0.38)
  for(i in 1:m) {
    for (j in 1:m) {
      if (i != j) {
        Sigma[i, j] <- rho^abs(i-j) * sqrt(Sigma[i, i] * Sigma[j, j])
      }
    }
  }
  random_effects <- rmvnorm(n = ncluster,  sigma = Sigma)
  
  data$u0 = rep(random_effects[, 1],
                each = 5, times=1)
  
  data$u1 = rep(random_effects[, 2],
                each = 5, times=1)
  
  sigmae = sqrt(0.45) # 0.6
  data$e = rnorm(N, 0, sigmae)
  
  data$b0 = -0.5
  data$b1 = 0.62
  
  data$Y8 = data$b0 + data$b1*data$Y5 + data$u1*data$Y5 + data$u0 + data$e
  
  
  data <- data[, !(names(data) %in% c("u0","u1", "e", "b0", "b1"))]
  
  #== Y9|Y3:Y4:Y5:G
  m=4
  rho <- 0.42
  Sigma <- matrix(NA, m, m)
  diag(Sigma) <- c(0.6, 0.52, 0.42, 0.35)
  for(i in 1:m) {
    for (j in 1:m) {
      if (i != j) {
        Sigma[i, j] <- rho^abs(i-j) * sqrt(Sigma[i, i] * Sigma[j, j])
      }
    }
  }
  random_effects <- rmvnorm(n = ncluster,  sigma = Sigma)
  
  data$u0 = rep(random_effects[, 1],
                each = 5, times=1)
  
  data$u1 = rep(random_effects[, 2],
                each = 5, times=1)
  data$u2 = rep(random_effects[, 3],
                each = 5, times=1)
  
  data$u3 = rep(random_effects[, 4],
                each = 5, times=1)
  sigmae = sqrt(0.55) # 0.85
  data$e = rnorm(N, 0, sigmae)
  
  data$b0 = 0.97
  data$b1 = 0.27
  data$b2 = 1.31
  data$b3 = -0.34
  
  data$Y9 = data$b0 + data$b1*data$Y3 +data$b2*data$Y4 +data$b3*data$Y5 + 
    data$u1*data$Y3 + data$u2*data$Y4 + data$u3*data$Y5 + 
    data$u0 + data$e
  

  data <- data[, !(names(data) %in% c("u0","u1","u2","u3", "e", "b0", "b1", "b2","b3"))]
  
  #== Y10|Y4:G
  m= 2
  rho <- 0.44
  Sigma <- matrix(NA, m, m)
  diag(Sigma) <- c(0.58, 0.42)
  for(i in 1:m) {
    for (j in 1:m) {
      if (i != j) {
        Sigma[i, j] <- rho^abs(i-j) * sqrt(Sigma[i, i] * Sigma[j, j])
      }
    }
  }
  random_effects <- rmvnorm(n = ncluster,  sigma = Sigma)
  
  data$u0 = rep(random_effects[, 1],
                each = 5, times=1)
  
  data$u1 = rep(random_effects[, 2],
                each = 5, times=1)
  
  sigmae = sqrt(0.52) # 0.6
  data$e = rnorm(N, 0, sigmae)
  
  data$b0 = 1.92
  data$b1 = 1.19
  
  data$Y10 = data$b0 + data$b1*data$Y4 + data$u1*data$Y4 + data$u0 + data$e
  
  data <- data[, !(names(data) %in% c("u0","u1", "e", "b0", "b1", "ID"))]
  
  
  # === True DAG and True Parameters
  #True DAG
  TMBN = model2network("[G][Y1|G][Y3|G][Y2|Y1:G][Y5|Y3:G][Y6|Y1:Y3:G][Y4|Y2:G][Y8|Y5:G][Y7|Y4:G][Y9|Y3:Y4:Y5:G][Y10|Y4:G]")
  
  #== Structure learning ==#
  data$G<- as.factor(data$G)
  inbound.arcs = data.frame(from = setdiff(names(data), "G"), to = "G")
  outbound.arcs = data.frame(from = "G", to = setdiff(names(data), "G"))
  
  
  dag.lme = hc(data, score = "custom-score", fun = scoreME, args = list(group = "G"),
               whitelist = outbound.arcs)
  dag.inla = hc(data, score = "custom-score", fun = scoreME1, args = list(group = "G"),
                whitelist = outbound.arcs)
  dag.inla1 = hc(data, score = "custom-score", fun = scoreME2, args = list(group = "G"),
                 whitelist = outbound.arcs)
  
  TMBN_cp <-cpdag(TMBN)
  dag.lme_cp<-cpdag(dag.lme)
  dag.inla_cp<-cpdag(dag.inla)
  dag.inla1_cp<-cpdag(dag.inla1)
  
  shd_lme[[k]] = shd(TMBN_cp, dag.lme_cp) 
  shd_inla[[k]] = shd(TMBN_cp, dag.inla_cp)
  shd_inla1[[k]] = shd(TMBN_cp, dag.inla1_cp) 
  
  
}

shd_lme = unlist(shd_lme) 
shd_inla = unlist(shd_inla) 
shd_inla1 = unlist(shd_inla1)  

stats_1<-data.frame(shd_lme, shd_inla, shd_inla1)

setwd("path")
out.file2 <- sprintf(
  "Gauss_MbyM2_%s_%s.rds",
  Sys.Date(), format(Sys.time(), "%H%M%S")
)
saveRDS(stats_1, file = out.file2)

