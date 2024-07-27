# OWN FUNS---------------------------
# Theil U RMSE
THEIL_U_RMSE <- function (x, y){
  a <- sqrt(sum((x)^2))
  b <- sqrt(sum((y)^2))
  output <- a/b
  return(output)
}


THEIL_U_RMSE_R <- function (x, y){
  a <- (sum((x)^2))
  b <- (sum((y)^2))
  output <- a/b
  return(output)
}


# -----------------------------------------------------------------------

# Theil U MAE

THEIL_U_MAE <- function (x, y){
  a <- sqrt(sum(abs(x)))
  b <- sqrt(sum(abs(y)))
  output <- a/b
  return(output)
}


# -----------------------------------------------------------------------
# Clark West test simple
# my implementation

clark_west <- function(actual, restricted, unrestricted,horizon){
  # restricted smaller model
  
  # f_hat <- e.m1^2-(e.m2^2-(yf.m1-yf.m2)^2)
  
  e1 = actual - restricted
  e2 = actual - unrestricted
  f_hat = e1^2 - (e2^2 - (restricted - unrestricted)^2)
  
  model_cw=lm(f_hat~1)
  
  beta_cw = coef(model_cw)[1]
  # getting estimate
  # standard  error
  # B_SE = sqrt(vcov(model_cw))
  
  # newey west error
  SEs <- NeweyWest(model_cw, lag = horizon, prewhite=FALSE)
  B_SE <- sqrt(diag(SEs))
  
  # extracting standard error
  CW=as.numeric(beta_cw/B_SE)
  pv =  1-pnorm(CW,0,1)
  
  
  
  results <- list(test=0,pvalue=0)
  results$test <- CW
  results$pvalue <- pv
  return(results)
}



# Clark West test
cw <- function(e.m1,e.m2,yf.m1,yf.m2,horizon){
  len=min(length(e.m1),length(e.m2),length(yf.m1),length(yf.m2))
  e.m1 <- e.m1[1:len]
  e.m2 <- e.m2[1:len]
  yf.m1 <- yf.m1[1:len]
  yf.m2 <- yf.m2[1:len]
  
  
  nw <- function(y,qn){
    #input: y is a T*k vector and qn is the truncation lag
    #output: the newey west HAC covariance estimator
    #Formulas are from Hayashi
    T <- length(y)
    ybar <- rep(1,T) * ((sum(y))/T)
    dy <- y-ybar
    G0 <- t(dy) %*% dy/T
    for (j in 1:qn){
      gamma <- t(dy[(j+1):T]) %*% dy[1:(T-j)]/(T-1)
      G0 <- G0+(gamma+t(gamma))*(1-abs(j/qn))
    }
    return(as.numeric(G0))
  }
  
  P <- length(e.m1)
  froll.adj <- e.m1^2-(e.m2^2-(yf.m1-yf.m2)^2)
  varfroll.adj <- nw(froll.adj,horizon)
  CW <- sqrt(P)*(mean(froll.adj))/sqrt(varfroll.adj)
  pv <- 1-pnorm(CW,0,1)
  results <- list(test=0,pvalue=0)
  results$test <- CW
  results$pvalue <- pv
  return(results)
}

# -----------------------------------------------------------------------
ewma <- function(x, v = 0.98, y0) {
  ## calculate exponentially-weighted moving average as in Cieslak-Povala
  ## v = 0.95 is the default smoothing parameter, appropriate for quarterly data
  tau <- numeric(length(x))
  if (missing(y0)) {
    tau[1] <- x[1]
  } else {
    tau[1] <- y0
  }
  for (t in 2:length(x))
    tau[t] <- tau[t-1] + (1-v)*(x[t] - tau[t-1])
  tau
}

covEWMA <- function(factors, lambda=0.96, return.cor=FALSE) {
  ## Inputs:
  ## factors    N x K numerical factors data.  data is class data.frame
  ##            N is the time length and K is the number of the factors.  
  ## lambda     scalar. exponetial decay factor between 0 and 1. 
  ## return.cor Logical, if TRUE then return EWMA correlation matrices
  ## Output:  
  ## cov.f.ewma  array. dimension is N x K x K.
  ## comments:
  ## 1. add optional argument cov.start to specify initial covariance matrix
  ## 2. allow data input to be data class to be any rectangular data object
  
  
  if (is.data.frame(factors)){
    factor.names  = colnames(factors)
    t.factor      = nrow(factors)
    k.factor      = ncol(factors)
    factors       = as.matrix(factors)
    t.names       = rownames(factors)
  } else {
    stop("factor data should be saved in data.frame class.") 
  }
  if (lambda>=1 || lambda <= 0){
    stop("exponential decay value lambda should be between 0 and 1.")
  } else {
    cov.f.ewma = array(,c(t.factor,k.factor,k.factor))
    cov.f = var(factors)  # unconditional variance as EWMA at time = 0 
    FF = (factors[1,]- mean(factors)) %*% t(factors[1,]- mean(factors))
    cov.f.ewma[1,,] = (1-lambda)*FF  + lambda*cov.f
    for (i in 2:t.factor) {
      FF = (factors[i,]- mean(factors)) %*% t(factors[i,]- mean(factors))
      cov.f.ewma[i,,] = (1-lambda)*FF  + lambda*cov.f.ewma[(i-1),,]
    }
    
  }
  # 9/15/11: add dimnames to array
  dimnames(cov.f.ewma) = list(t.names, factor.names, factor.names)
  
  if(return.cor) {
    cor.f.ewma = cov.f.ewma
    for (i in 1:dim(cor.f.ewma)[1]) {
      cor.f.ewma[i, , ] = cov2cor(cov.f.ewma[i, ,])
    }
    return(cor.f.ewma)
  } else{
    return(cov.f.ewma)  
  }
}


# Bootstrap p value
# Maio (2013)
pval_boot_fun <- function(x=monthly$lty_detr, y=monthly$equity_premium_lead,replication=10000) {
  
  y <- as.numeric(y)
  x <- as.numeric(x)
  
  y <- y[-1]
  reg_res1 <- (lm(y~ x[-1]))
  reg_res2 <- (lm(x~ lag(x)))
  
  summary(reg_res1)
  
  residuals1 <- residuals(reg_res1)
  residuals2 <- residuals(reg_res2)
  
  
  coef_sim <- c()
  
  for (i in 1:replication) {
    kontrola=NA
    while(is.na(kontrola)) {  
      sample_round <- as.numeric(sample(1:length(residuals2), replace = TRUE))
      x_sample <- as.numeric(sample(residuals2,1))
      x_boot <- x_sample
      y_boot <- c()
      for (j in 1:length(y)) {
        x_sample <- as.numeric(coef(reg_res2)[1]+ coef(reg_res2)[2]*x_sample +residuals2[sample_round[j]])
        y_sample <-  as.numeric(coef(reg_res1)[1]+residuals1[sample_round[j]])
        
        x_boot <- rbind(x_boot,x_sample)
        y_boot <- rbind(y_boot,y_sample)
        
      }
      x_boot <- x_boot[-1]
      reg_res_sim <- (lm(y_boot~ x_boot))
      
      kontrola <- coef(reg_res_sim)[2]
    }
    coef_sim <- rbind(coef_sim,coef(reg_res_sim)[2])
    
    
    
  }
  
  
  
  
  
  if(as.numeric(coef(reg_res1)[2])<0){
    res <- table((as.numeric(coef_sim)<=as.numeric(coef(reg_res1)[2])))[2]+
      table((as.numeric(coef_sim)>=-as.numeric(coef(reg_res1)[2])))[2]
    
  }
  if(as.numeric(coef(reg_res1)[2])>0){
    res <- table((as.numeric(coef_sim)>=as.numeric(coef(reg_res1)[2])))[2]+
      table((as.numeric(coef_sim)<=-as.numeric(coef(reg_res1)[2])))[2]
    
  }
  
  
  pval_boot <- res/replication
  return(pval_boot)
}

# ---------------------
