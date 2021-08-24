#Estimadores insesgados

first <- function(clusters, Mi) {
  m <- c()
  tau <- c()
  mu <- c()
  S <- c()
  for (i in 1:length(clusters)) {
    m[i] <- length(clusters[[i]])
    tau[i] <- Mi[i]/m[i]*sum(clusters[[i]])
    mu[i] <- tau[i]/Mi[i]
    S[i] <- var(clusters[[i]])
  }
  return(list(mi = m, taui = tau, mui = mu, Si2 = S))
}

tau_mu_2 <- function(taui, Si2, Mi, mi, N) {
  n <- length(taui)
  Mo <- sum(Mi)
  tau2 <- (N/n)*sum(taui)
  mu2 <- tau2/Mo
  taubar <- mean(taui)
  supm2 <- sum((taui - taubar)^2)/(n-1)
  vartau2 <- N^2*(1-n/N)*supm2/n + (N/n)*sum(Mi^2*(1-mi/Mi)*Si2/mi)
  varmu2 <- vartau2/Mo^2
  Btau2 <- 2*sqrt(vartau2)
  Bmu2 <- 2*sqrt(varmu2)
  Litau2 <- tau2 - Btau2; Lstau2 <- tau2 + Btau2
  Limu2 <- mu2 - Bmu2; Lsmu2 <- mu2 + Bmu2
  result <- data.frame(estimations = c(tau2, mu2), B = c(Btau2, Bmu2),
                       Li = c(Litau2, Limu2), Ls = c(Lstau2, Lsmu2))
  return(result)
}





#Estimadores de razón



#PPT



#Conglomerados igual tamaño y tamaño de muestra