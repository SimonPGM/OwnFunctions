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

tau_mu_2 <- function(taui, Si2, Mi, mi, N, Mo) {
  n <- length(taui)
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
  rownames(result) <- c("tau(2)", "mu(2)")
  return(result)
}

#Estimadores de razón

EstRazMu <- function(Mi, N, n, mi, si2,  ti = NULL, yibar = NULL, M0 = NULL){
  #Mi Tamaños de los conglomerados
  #N Número total de conglomerados
  #n Número de conglomerados en la muestra
  #mi Tamaño de muestra en el i-ésimo conglomerado
  #si2 <- Varianza dentro de los conglomerados
  #Deben dar ti o yibar (totales o medias por conglomerado) 
  if (is.null(M0)) {
    M0 <- N*mean(Mi)
  }
  Mbar <- M0/N
  if (is.null(ti)){
    Mu_r2 <- sum(Mi*yibar)/sum(Mi)
    s2_rupm <- sum(Mi^2*(yibar-Mu_r2)^2/(n-1))
  }
  else{
    Mu_r2 <- sum(ti)/sum(Mi)
    s2_rupm <- sum((ti-Mi*Mu_r2)^2/(n-1))
  }
  varMu_r2 <- ((N-n)/N)*(1/Mbar^2)*(s2_rupm/n)+(1/(n*N*Mbar^2))*sum(Mi^2*((Mi-mi)/Mi)*(si2/mi))
  t_r2 <- M0*Mu_r2
  vart_r2 <- M0^2*varMu_r2
  LEE <- c(2*sqrt(varMu_r2), 2*sqrt(vart_r2)) 
  intervalo <- data.frame(LI = c(Mu_r2 - LEE[1], t_r2 - LEE[2]), LS = c(Mu_r2 + LEE[1], t_r2 + LEE[2]))
  rownames(intervalo) <- c("Mu_r2", "T_r2")
  list(Estimaciones = data.frame(Mu_r2 = Mu_r2, T_r2 = t_r2), 
       Varianzas = data.frame(varMu_r2 = varMu_r2, varT_r2 = vart_r2), 
       LEEs = data.frame(LEEMu_r2 = LEE[1], LEET_r2 = LEE[2]), 
       Intervalo = intervalo, S2_rUPM = s2_rupm)
}

EstRazProp <- function(Mi, N, n, mi, pi, M0 = NULL){
  #Mi Tamaños de los conglomerados
  #N Número total de conglomerados
  #n Número de conglomerados en la muestra
  #mi Tamaño de muestra en el i-ésimo conglomerado
  if (is.null(M0)) {
    M0 <- N*mean(Mi)
  }
  Mbar <- M0/N
  p_r2 <- sum(Mi*pi)/sum(Mi)
  s2_rupm <- sum(Mi^2*(pi-p_r2)^2/(n-1))
  varp_r2 <- ((N-n)/N)*(1/Mbar^2)*(s2_rupm/n)+(1/(n*N*Mbar^2))*sum(Mi^2*((Mi-mi)/Mi)*(pi*(1-pi)/(mi-1)))
  A_r2 <- M0*p_r2
  varA_r2 <- M0^2*varp_r2
  LEE <- c(2*sqrt(varp_r2), 2*sqrt(varA_r2)) 
  intervalo <- data.frame(LI = c(p_r2 - LEE[1], A_r2 - LEE[2]), LS = c(p_r2 + LEE[1], A_r2 + LEE[2]))
  rownames(intervalo) <- c("P_r2", "A_r2")
  list(Estimaciones = data.frame(P_r2 = p_r2, A_r2 = A_r2), 
       Varianzas = data.frame(varP_r2 = varp_r2, varA_r2 = varA_r2), 
       LEEs = data.frame(LEEP_r2 = LEE[1], LEEA = LEE[2]), 
       Intervalo = intervalo, S2_rUPM = s2_rupm)
}

#PPT

Estimacionesppt2 <- function(Mi, n, M0, ti = NULL, yibar = NULL){
  #Mi Tamaños de los conglomerados
  #n Número de conglomerados en la muestra
  #deben dar ti o yibar (totales o medias por conglomerado)
  y_ibar <-  ifelse(rep(is.null(ti), n), yibar, ti/Mi)
  
  
  t_ppt2 <- M0/n*sum(y_ibar)
  mu_ppt2 <- t_ppt2/M0
  vart_ppt2 <- M0^2/(n*(n-1))*sum((y_ibar-mu_ppt2)^2)
  varmu_ppt2 <- vart_ppt2/M0^2
  intervalos <- data.frame(LI = c(t_ppt2 - 2*sqrt(vart_ppt2), mu_ppt2 - 2*sqrt(varmu_ppt2)),
                           LS = c(t_ppt2 + 2*sqrt(vart_ppt2), mu_ppt2 + 2*sqrt(varmu_ppt2)))
  rownames(intervalos) <- c("T_ppt2", "Mu_ppt2")
  LEE <- data.frame(T_ppt2 = 2*sqrt(vart_ppt2), Mu_ppt2 = 2*sqrt(varmu_ppt2))
  vars <- data.frame(T_ppt2 = vart_ppt2, Mu_ppt2 = varmu_ppt2)
  list(Estimaciones = data.frame(T_ppt2 = t_ppt2, Mu_ppt2 = mu_ppt2),
       Mu_i = y_ibar, VAR = vars, LEE = LEE, Intervalos = intervalos)
}

#Polietápico con Reemplazo

EstPoliReem <- function(Mi, n, M0, ti = NULL,  yibar = NULL){
  #Mi Tamaños de los conglomerados
  #n Número de conglomerados en la muestra
  #deben dar ti o yibar (totales o medias por conglomerado)
  pi <- Mi/M0
  if (is.null(yibar)){
    t_p <- M0/n*sum(ti/Mi)
    vart_p <- 1/(n*(n-1))*sum((ti/pi-t_p)^2)
  }
  else{
  t_p <- M0/n*sum(yibar) 
  vart_p <- 1/(n*(n-1))*sum(((Mi*yibar)/pi-t_p)^2)}
  data.frame(T_p = t_p, Var = vart_p)
}

#Conglomerados igual tamaño y tamaño de muestra
biphase_equal_size_mu <- function(yi_bar, Si2, N, n, M, m){
  #Estimacion puntual
  mu2 <- 1/n * sum(yi_bar)
  
  #Varianza
  #Si2 son varianzas 
  #mu2 = mu_hat
  MSB <- m/(n - 1) * sum((yi_bar - mu2)^2)
  MSW <- 1/n * sum(Si2)
  var_mu2 <- (1 - n/N) * MSB/(n*m) + (1 - m/M) * 1/N * MSW/m
  
  #Error estandar
  ee_mu2 <- sqrt(var_mu2)
  
  #Limite de error de estimacion
  B_mu2 <- 2 * ee_mu2
  
  #IC
  IC <- c(mu2 - B_mu2, mu2 + B_mu2)
  
  overall <- list(mu_2_hat = mu_2,
                  Variance = var_mu2,
                  Standard_Error = ee_mu2,
                  LEE_B = B_mu2,
                  IC = IC)
  return(overall)
}


sample_size_biphase_equal <- function(yi_bar, Si2, n_pilot, m_pilot, 
                                      ci, c = NULL, var_mu2 = NULL){
  #Se necesita valores inciales de n, m Si2 de una muestra piloto
  mu2 <- 1/n_pilot * sum(yi_bar)
  MSB <- m_pilot/(n_pilot - 1) * sum((yi_bar - mu2)^2)
  MSW <- 1/n_pilot * sum(Si2) #Recordar que MSW estima a sigma2w_hat
  sigma2b_hat <- 1/m_pilot *(MSB - MSW)
  f1 <- MSW/sigma2b_hat
  f2 <- ci[1]/ci[2]
  m <- sqrt(f1 * f2)
  if(!is.null(c)){
    n <- c/(ci[1] + round(m) * ci[2])
  }
  else{
    n <- 1/var_mu2 * (sigma2b_hat + MSW/m)
  }
  overall <- data.frame(actual = c(m, n),
                        approx = round(c(m, n)))
  rownames(overall) <- c("m", "n")
  return(overall)
}




