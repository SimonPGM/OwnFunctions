######Normalidad
mu_nosigma <- function(mu0, datos, alpha){
  #Function for ht in a multivariate sample with normal distribution when varcov matrix is unknown
  mu0 <- matrix(mu0, ncol = 1)
  xbar <- apply(datos, 2, mean)
  xbar <- matrix(xbar, ncol = 1)
  n <- dim(datos)[1]
  p <- dim(datos)[2]
  S <- var(datos)
  k <- (n-1)*p/(n-p)
  f <- (n/k)*t(xbar-mu0) %*% solve(S) %*% (xbar-mu0)
  pvalue <- pf(f, df1 = p, df2 = n-p, lower.tail = F)
  reject <- f > qf(alpha, df1 = p, df2 = n-p, lower.tail = F)
  list(pvalue = pvalue, f0 = f, reject = reject)
}

mu_sigma <- function(mu0, datos, sigma, alpha) {
  #Function for ht in a multivariate sample with normal distribution when varcov matrix is known
  xbar <- as.matrix(apply(datos, 2, mean), ncol = 1)
  mu0 <- as.matrix(mu0, ncol = 1)
  p <- dim(datos)[2]
  chi0 <- n*t(xbar-mu0) %*% solve(sigma) %*% (xbar - mu0)
  pvalue <- pchisq(chi0, df = p, lower.tail = F)
  reject <- chi0 > qchisq(alpha, df = p, lower.tail = F)
  list(pvalue = pvalue, chi0 = chi0, reject = reject)
}

two_mu_sigmaequals_unknown <- function(datos1, datos2, alpha, delta0 = NULL) {
  #difference of means equals to delta0 if varcovs are unknow and equals
  p <- dim(datos1)[2]
  n <- dim(datos1)[1]
  m <- dim(datos2)[1]
  
  if (is.null(delta0)) {
    delta0 <- matrix(rep(0, p), ncol = 1)
  } else {
    delta0 <- matrix(delta0, ncol = 1)
  }
  
  xbar <- matrix(apply(datos1, 2, mean), ncol = 1)
  ybar <- matrix(apply(datos2, 2, mean), ncol = 1)
  S1 <- var(datos1)
  S2 <- var(datos2)
  Sp <- ((n-1)*S1 + (m-1)*S2)/(n+m-2)
  k <- (n+m-2)*p/(n+m-p-1)
  f0 <- (n*m)/((n+m)*k)*t(xbar - ybar - delta0) %*% solve(Sp) %*% (xbar - ybar - delta0)
  pvalue <- pf(f0, p, n+m-p-1, lower.tail = F)
  reject <- f0 > qf(alpha, p, n+m-p-1, lower.tail = F)
  list(pvalue = pvalue, f0 = f0, reject = reject)
  
}

two_mu_sigmanotequals_unknown <- function(datos1, datos2, alpha, delta0 = NULL) {
  
  p <- dim(datos1)[2]
  n <- dim(datos1)[1]
  m <- dim(datos2)[1]
  
  f <- function(X) {
    n <- dim(X)[1]
    (sum(diag(X)) + sum(diag(X))^2)/(n-1) 
  }
  
  if(is.null(delta0)) {
    delta0 <- matrix(rep(0, p), ncol = 1)
  } else {
    delta0 <- matrix(delta0, ncol = 1)
  }
  
  xbar <- matrix(apply(datos1, 2, mean), ncol = 1)
  ybar <- matrix(apply(datos2, 2, mean), ncol = 1)
  V1 <- var(datos1)/n
  V2 <- var(datos2)/m
  Se <- V1 + V2
  num <- (sum(diag(Se)) + sum(diag(Se))^2)
  den <- sum(unlist(lapply(list(V1, V2), f)))
  v <- round(num/den)
  k <- v*p/(v-p+1)
  f0 <- (1/k)*t(xbar-ybar-delta0) %*% solve(Se) %*% (xbar-ybar-delta0)
  pvalue <- pf(f0, p, v-p+1, lower.tail = F)
  reject <- f0 > qf(alpha, p, v-p+1, lower.tail = F)
  list(pvalue = pvalue, f0 = f0, reject = reject)
  
}

two_mu_sigmaequals_known <- function(datos1, datos2, sigma, alpha, delta0 = NULL) {
  
  p <- dim(datos1)[2]
  n <- dim(datos1)[1]
  m <- dim(datos2)[1]
  
  if(is.null(delta0)) {
    delta0 <- matrix(rep(0, p), ncol = 1)
  } else {
    delta0 <- matrix(delta0, ncol = 1)
  }

  xbar <- matrix(apply(datos1, 2, mean), ncol = 1)
  ybar <- matrix(apply(datos2, 2, mean), ncol = 1)
  
  chi0 <- n*m/(n+m) * t(xbar-ybar-delta0) %*% solve(sigma) %*% (xbar-ybar-delta0)
  pvalue <- pchisq(chi0, p, lower.tail = F)
  reject <- chi0 > qchisq(alpha, p, lower.tail = F)
}

two_mu_paired_unknownsigma <- function(datos1, datos2, sigma, alpha) {
  
  D <- datos1 - datos2
  Dbar <- matrix(apply(D, 2, mean), ncol = 1)
  SD <- var(D)
  p <- dim(D)[2]
  n <- dim(D)[1]
  
  k <- (n-1)*p/(n-p)
  f0 <- (n/k)*t(Dbar) %*% solve(S) %*% Dbar
  pvalue <- pf(f0, p, n-p, lower.tail = F)
  reject <- f0 > qf(alpha, p, n-p, lower.tail = F)
  list(pvalue = pvalue, f0 = f0, reject = reject)
  
}

contrast_sigma_unknown <- function(datos, C, gamma, alpha) {
  
  k <- length(gamma)
  n <- dim(datos)[1]
  
  gamma <- matrix(gamma, ncol = 1)
  xbar <- matrix(apply(datos, 2, mean), ncol = 1)
  unb <- C %*% xbar
  S <- var(datos)
  varcov <- C %*% S %*% t(C)
  
  c <- (n-1)*k/(n-k)
  
  f0 <- (n/c)*t(unb-gamma) %*% solve(varcov) %*% (unb-gamma)
  pvalue <- pf(f0, k, n-k, lower.tail = F)
  reject <- f0 > qf(alpha, k, n-k, lower.tail = F)
  list(pvalue = pvalue, f0 = f0, reject = reject)
  
  
}

contrast_sigma_known <- function(datos, C, gamma, sigma, singif) {
  
  k <- length(gamma)
  n <- dim(datos)[1]
  
  gamma <- matrix(gamma, ncol = 1)
  xbar <- matrix(apply(datos, 2, mean), ncol = 1)
  unb <- C %*% xbar
  varcov <- C %*% sigma  %*% t(C)
  
  chi0 <- n*t(unb-gamma) %*% solve(varcov) %*% (unb-gamma)
  pvalue <- pchisq(chi0, k, lower.tail = F)
  reject <- chi0 > qchisq(alpha, k, lower.tail = F)
  list(pvalue = pvalue, chi0 = chi0, reject = reject)
  
}

contrast_tlc <- function(datos, C, gamma, alpha) {
  
  m <- dim(C)[1]
  n <- dim(datos)[1]
  
  gamma <- matrix(gamma, ncol = 1)
  xbar <- matrix(apply(datos, 2, mean), ncol = 1)
  unb <- C %*% xbar
  S <- var(datos)
  varcov <- C %*% S %*% t(C)
  
  
  chi0 <- n*t(unb-gamma) %*% solve(varcov) %*% (unb-gamma)
  pvalue <- pchisq(chi0, m, lower.tail = F)
  reject <- chi0 > qchisq(alpha, m, lower.tail = F)
  list(pvalue = pvalue, chi0 = chi0, reject = reject)
  
}

M_test <- function(listadatos, alpha) {
  
  g <- length(listadatos)
  p <- dim(listadatos[[1]])[2]
  
  Sl <- lapply(listadatos, var)
  vl <- unlist(lapply(listadatos, function(x) dim(x)[1] -1 ))
  v <- sum(vl)
  temp1 <- 0
  temp2 <- 0
  
  for (i in 1:g) {
    temp1 <- vl[i]*Sl[[i]] + temp1
    temp2 <- vl[i]*log(det(Sl[[i]])) + temp2
  }
  
  Sp <- (1/v)*temp1
  M <- v*log(det(Sp)) - temp2
  
  u <- (sum(1/vl) - 1/v)*(2*p^2 + 3*p - 1 )/(6*(p+1)*(g-1))
  C <- (1-u)*M
  k <- (g-1)*p*(p+1)/2
  pvalor <- pchisq(C, k, lower.tail = F)
  reject <- C > qchisq(alpha, k, lower.tail = F)
  list(pvalor = pvalor, C = C, reject = reject)
}


#CENTRAL LIMIT THEOREM AND LRT
sigma_normal_lrt <- function(sigma0, datos, alpha, aprox = T){
  n <- dim(datos)[1]; p <- dim(datos)[2]
  S <- var(datos)
  tr <- sum(diag(S %*% solve(sigma0)))
  if(!aprox){
    aux1 <- ((1 - 1/n)^p * (det(S)/det(sigma0)))^(n/2) 
    aux2 <- exp(-1/2 * ((n - 1) * tr - n*p))
  }
  else{
    aux1 <- (det(S)/det(sigma0))^(n/2)
    aux2 <- exp(-n/2 * (tr - p))
  }
  lambda <- aux1 * aux2
  lambda_star <- -2*log(lambda)
  k <- p*(p + 1)/2
  if(n > 30){
    RR <- lambda_star > qchisq(alpha, k, lower.tail = F)
    pvalue <- pchisq(lambda_star, k, lower.tail = F)
  }
  else{
    correct <- 1 - 1/(6*(n - 1)) * (2*p + 1 - 2/(p + 1))
    lambda_star <- correct * lambda_star
    RR <- lambda_star > qchisq(alpha, k, lower.tail = F) 
    pvalue <- pchisq(lambda_star, k, lower.tail = F)
  }
  overall <- list(pvalue = pvalue,
                  Lambda_star = lambda_star,
                  Quantile = qchisq(alpha, k, lower.tail = F),
                  Reject = RR)
  return(overall)
}

#PH para mu usando TLC
mu_nosigma_lct <- function(mu0, datos, alpha){
  if(!is.matrix(mu0)){
    mu0 <- matrix(mu0, ncol = 1)
  }
  n <- dim(datos)[1]; p <- dim(datos)[2]
  x_bar <- apply(datos, 2, mean) %>% as.matrix(ncol = 1)
  S <- var(datos)
  Chi_02 <- n * t(x_bar - mu0) %*% solve(S) %*% (x_bar - mu0)
  RR <- Chi_02 > qchisq(alpha, p, lower.tail = F)
  pvalue <- pchisq(Chi_02, p, lower.tail = F)
  overall <- list(pvalue = pvalue, Chi_Squared_Statistic = Chi_02,
                  Quantile = qchisq(alpha, p, lower.tail = F),
                  Reject = RR)
  return(overall)
} 


#FUNCIONES PARA DIFERENCIAS DE MEDIAS USANDO TLC
two_mu_lct_sigmaequals_unknown <- function(datos_x, datos_y, alpha, gamma0 = NULL){
  n <- dim(datos_x)[1]; m <- dim(datos_y)[1]; p <- dim(datos_x)[2]
  x_bar <- apply(datos_x, 2, mean) %>% as.matrix(ncol = 1)
  y_bar <- apply(datos_y, 2, mean) %>% as.matrix(ncol = 1)
  Sx <- var(datos_x); Sy <- var(datos_y)
  Sp <- ((n - 1)*Sx + (m - 1)*Sy)/(n + m - 2)
  if(is.null(gamma0)){
    Chi_02 <- (n*m/(n + m)) * t(x_bar - y_bar) %*% solve(Sp) %*% (x_bar - y_bar)
  }
  else{
    if(!is.matrix(gamma0)){
      gamma0 <- matrix(gamma0, ncol = 1)
    }
    Chi_02 <- (n*m/(n + m)) * t(x_bar - y_bar - gamma0) %*% solve(Sp) %*% 
      (x_bar - y_bar - gamma0)
  }
  RR <- Chi_02 > qchisq(alpha, p, lower.tail = F)
  pvalue <- pchisq(Chi_02, p, lower.tail = F)
  overall <- list(pvalue = pvalue, Chi_Squared_Statistic = Chi_02,
                  Quantile = qchisq(alpha, p, lower.tail = F),
                  Reject = RR)
  return(overall)
}


two_mu_lct_sigmanotequals_unknown <- function(datos_x, datos_y, alpha, gamma0 = NULL){
  n <- dim(datos_x)[1]; m <- dim(datos_y)[1]; p <- dim(datos_x)[2]
  x_bar <- apply(datos_x, 2, mean) %>% as.matrix(ncol = 1)
  y_bar <- apply(datos_y, 2, mean) %>% as.matrix(ncol = 1)
  Sx <- var(datos_x); Sy <- var(datos_y)
  Se <- Sx/n + Sy/m
  if(is.null(gamma0)){
    Chi_02 <- t(x_bar - y_bar) %*% solve(Se) %*% (x_bar - y_bar)
  }
  else{
    if(!is.matrix(gamma0)){
      gamma0 <- matrix(gamma0, ncol = 1)
    }
    Chi_02 <- t(x_bar - y_bar - gamma0) %*% solve(Se) %*% (x_bar - y_bar - gamma0)
  }
  RR <- Chi_02 > qchisq(alpha, p, lower.tail = F)
  pvalue <- pchisq(Chi_02, p, lower.tail = F)
  overall <- list(pvalue = pvalue, Chi_Squared_Statistic = Chi_02,
                  Quantile = qchisq(alpha, p, lower.tail = F),
                  Reject = RR)
  return(overall)
}


two_mu_lct_sigmaequals_known <- function(datos_x, datos_y, Sigma, alpha, gamma0 = NULL){
  n <- dim(datos_x)[1]; m <- dim(datos_y)[1]; p <- dim(datos_x)[2]
  x_bar <- apply(datos_x, 2, mean) %>% as.matrix(ncol = 1)
  y_bar <- apply(datos_y, 2, mean) %>% as.matrix(ncol = 1)
  if(is.null(gamma0)){
    Chi_02 <- (n*m/(n + m)) * t(x_bar - y_bar) %*% solve(Sigma) %*% (x_bar - y_bar)
  }
  else{
    if(!is.matrix(gamma0)){
      gamma0 <- matrix(gamma0, ncol = 1)
    }
    Chi_02 <- (n*m/(n + m)) * t(x_bar - y_bar - gamma0) %*% solve(Sigma) %*% 
      (x_bar - y_bar - gamma0)
  }
  RR <- Chi_02 > qchisq(alpha, p, lower.tail = F)
  pvalue <- pchisq(Chi_02, p, lower.tail = F)
  overall <- list(pvalue = pvalue, Chi_Squared_Statistic = Chi_02,
                  Quantile = qchisq(alpha, p, lower.tail = F),
                  Reject = RR)
  return(overall)
}

#Regiones e intervalos simult?neos de confianza

reg_conf_mu <- function(mu0, means, s, n, alpha){
  #mu0 es el vector de H0
  #means es el vector de medias muestral
  #s es la matriz de covarianzas muestral
  #n es n xd 
  #alpha es la significancia del hiperelipsoide 
  
  p <- length(mu0)
  mu0 <- matrix(mu0, byrow = T, nrow = p)
  means <- matrix(means, byrow = T, nrow = p)
  t02 <- n*t(means-mu0)%*%solve(s)%*%(means-mu0)
  c2 <- (p*(n-1)/(n-p))*qf(alpha, p, n-p, lower.tail = F)
  eig <- eigen(s)
  
  SemiL <- c()
  for (i in 1:(p)) {
    SemiL[i] <- sqrt(eig$values[i])*1/sqrt(n)*sqrt(c2)
  }
  
  list(C2= c2, T02= t02,
       Resultado= ifelse(t02 > c2,"Mu0 no est? en la regi?n de confianza, se rechaza H0", 
                         "Mu0 est? en la regi?n de confianza, no se rechaza H0"),
       CentroElipsoide= means, SemilongitudEjes = SemiL, Ejes = eig$vectors)
}

int_sim_mu <- function(means, s, n, alpha){
  #means es el vector de medias muestral
  #s es la matriz de covarianza muestral
  #n es n xd x2
  #alpha es el nivel de confianza simult?neo
  
  p <- length(means)
  lower <- c()
  upper <- c()
  for (i in 1:p) {
    lower[i] <- means[i] - sqrt((p*(n-1)/(n-p))*qf(alpha, p, n-p, lower.tail = F))*sqrt(s[i,i]/n)
    upper[i] <- means[i] + sqrt((p*(n-1)/(n-p))*qf(alpha, p, n-p, lower.tail = F))*sqrt(s[i,i]/n)
  }
  
  cbind("L?mite Inferior" = lower, "L?mite Superior" = upper)
}

int_sim_diffmu <- function(means, s, n, index, alpha){
  #means es el vector de medias muestral
  #s es la matriz de covarianza muestral
  #n es n xd x3
  #index es un vector con los ?ndices de las diferencias 
  #alpha es el nivel de confianza simult?neo
  
  p <- length(means)
  i1 <- index[1]
  i2 <-index[2]
  low <- (means[i1]-means[i2])-sqrt((p*(n-1)/(n-p))*qf(alpha, p, n-p, lower.tail = F))*sqrt((s[i1,i1]+s[i2,i2]-2*s[i1,i2])/n) 
  up <- (means[i1]-means[i2])+sqrt((p*(n-1)/(n-p))*qf(alpha, p, n-p, lower.tail = F))*sqrt((s[i1,i1]+s[i2,i2]-2*s[i1,i2])/n) 
  cbind("L?mite Inferior" = low, "L?mite Superior" = up)
}

int_sim_bonf <- function(means, s, n, alpha){
  #means es el vector de medias muestral
  #s es la matriz de covarianza muestral
  #n es n xd x4
  #alpha es el nivel de confianza simult?neo
  
  p <- length(means)
  lower <- c()
  upper <- c()
  for (i in 1:p){
    lower[i] <- means[i] + qt(alpha/(2*p),n-1)*sqrt(s[i,i]/n)
    upper[i] <- means[i] - qt(alpha/(2*p),n-1)*sqrt(s[i,i]/n)
  }
  
  cbind("L?mite Inferior" = lower, "L?mite Superior" = upper)
}