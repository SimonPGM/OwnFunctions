#Conglomerados con tamaños iguales
mu_tau_con <- function(tau_i, M, N, n){
  #El tau_i, M y n se sacan de la BD
  #Si approx es FALSE de be especificar el nivel de confianza
  #Estimaciones puntuales
  mu_c_hat <- 1/n * sum(tau_i)
  mu_hat <- mu_c_hat/M
  tau_hat <- N * mu_c_hat
  
  #S2_con
  S2_con <- (1/(n - 1)) * sum((tau_i - mu_c_hat)^2) 
  
  #Varianzas estimadas de los estimadores
  var_mu_c_hat <- (1 - n/N) * (S2_con/n) 
  var_mu_hat <- (1/M^2) * var_mu_c_hat
  var_tau_hat <- N^2 * var_mu_c_hat
  Variances <- data.frame(var_mu_c_hat = var_mu_c_hat,
                          var_mu_hat = var_mu_hat,
                          var_tau_hat = var_tau_hat)
  #Errores estandar
  EE <- data.frame(ee_mu_c_hat = sqrt(var_mu_c_hat),
                   ee_mu_hat = sqrt(var_mu_hat),
                   ee_tau_hat = sqrt(var_tau_hat))
  
  #Limites de error de estimacion
  B_mu_c <- 2 * sqrt(var_mu_c_hat)
  B_mu_hat <- 2 * sqrt(var_mu_hat)
  B_tau_hat <- 2 * sqrt(var_tau_hat)
  LEE <- data.frame(B_mu_c = B_mu_c,
                    B_mu_hat = B_mu_hat,
                    B_tau_hat = B_tau_hat)
  
  #Intervalos de confianza
  IC_mu_c <- c(mu_c_hat - B_mu_c, mu_c_hat + B_mu_c)
  IC_mu_hat <- c(mu_hat - B_mu_hat, mu_hat + B_mu_hat)
  IC_tau_hat <- c(tau_hat - B_tau_hat, tau_hat + B_tau_hat)
  aux <- rbind(IC_mu_c, IC_mu_hat, IC_tau_hat)
  IC <- data.frame(lower = aux[, 1], upper = aux[, 2])
  rownames(IC) <- c("mu_c", "mu_hat", "tau_hat")
  
  #Resumen
  overall <- list(mu_c_hat = mu_c_hat, 
                  mu_hat = mu_hat,
                  tau_hat = tau_hat,
                  Variance = Variances,
                  Standard_Error = EE, 
                  LEE = LEE,
                  ICs = IC,
                  S2con = S2_con)
  return(overall)
}


con_equals_sample_size <- function(N, S2, D){
  #Recordar que el D varia dependiendo si el tamaño de muestra
  #es para mu o tau
  #mu: D = (B*M/Z)^2
  #tau: D = (B/(N*Z))^2
  num <- N*S2
  den <- N*D + S2
  n <- num/den
  overall <- data.frame(n_exacto = n, 
                        n_aproximado = ceiling(n))
  return(overall)
}


prop_tot_con <- function(M, N, n, Pi = NULL, Ai = NULL){
  #Si Pi es nulo, Ai NO puede ser nulo y viceversa
  #Solo es necesario ingresar alguna de las dos (Pi o Ai)
  if(is.null(Pi)){
    Pi <- Ai/M
  }
  #Estimacion puntual
  p_con <- 1/n * sum(Pi)
  A_con <- M*N*p_con
  
  #Varianzas estimadas
  var_p_con <- (1 - n/N) * (1/n) * (1/(n - 1)) * sum((Pi - p_con)^2)
  var_A_con <- (M*N)^2 * var_p_con
  varianzas <- data.frame(var_p_con = var_p_con,
                          var_A_con = var_A_con)
  
  #Errores estandars
  EE <- data.frame(ee_B_p = sqrt(var_p_con),
                   ee_B_A = sqrt(var_A_con))
  
  #Limites de error de estimacion
  B_p <- 2 * sqrt(var_p_con)
  B_A <- 2 * sqrt(var_A_con)
  LEE <- data.frame(B_p = B_p, B_A = B_A)
  
  #Intervalos de confianza
  IC_p <- c(p_con - B_p, p_con + B_p)
  #Para el total, el limite inferior se redondea por debajo
  #y el superior por encima (piso y techo respectivamente)
  IC_A <- c(floor(A_con - B_A), ceiling(A_con + B_A))
  aux <- rbind(IC_p, IC_A)
  IC <- data.frame(lower = aux[, 1],
                   upper = aux[, 2])
  rownames(IC) <- c("Prop", "Tot")
  
  #Resumen
  overall <- list(p_con = p_con, A_con = A_con,
                  Variance = varianzas,
                  Standard_Error = EE,
                  LEE = LEE,
                  ICs = IC)
  return(overall)
}

#Conglomerados tamaños distintos MAS y Razon
tau_muc_mu <- function(Mi, taui, N, Mo) {
  
  #Mi tamaño de los conglomerados (seleccionados en la muestra)
  #taui totales de los conglomerados (seleccionados en la muestra)
  #N total de conglomerados
  #Mo es el total de unidades muestrales en la poblacion
  
  n <- length(Mi)
  tauc.hat <- (N/n)*sum(taui)
  muc.hat <- tauc.hat/N
  M.bar <- Mo/N
  mu.hat <- muc.hat/M.bar
  temp <- (taui-muc.hat)^2
  vartauc.hat <- N^2*(1-n/N)*sum(temp)/(n*(n-1))
  varmuc.hat <- vartauc.hat/N^2
  varmu.hat <- vartauc.hat/Mo^2
  Btauc.hat <- 2*sqrt(vartauc.hat)
  Bmuc.hat <- 2*sqrt(varmuc.hat)
  Bmu.hat <- 2*sqrt(varmu.hat)
  Litauc <- tauc.hat - Btauc.hat; Lstauc <- tauc.hat + Btauc.hat
  Limuc <- muc.hat - Bmuc.hat; Lsmuc <- muc.hat + Bmuc.hat
  Limu <- mu.hat - Bmu.hat; Lsmu <- mu.hat + Bmu.hat
  result <- data.frame(Estimation = c(tauc.hat, muc.hat, mu.hat),
                       B = c(Btauc.hat, Bmuc.hat, Bmu.hat),
                       LI = c(Litauc, Limuc, Limu),
                       LS = c(Lstauc, Lsmuc, Lsmu))
  rownames(result) <- c("Tau", "Mu_c", "Mu")
  S2con <- sum(temp)/(n-1)
  return(list(result = result, S2con = S2con))
}

tau_muc_mu_r <- function(Mi, taui, N, Mo = NULL) {
  
  #Mi tamaño de los conglomerados
  #Taui totales de los conglomerados
  #N total de conglomerados
  #n número de coglomerados seleccionados en la muestra
  #Mo es el total de unidades muestrales en la poblacion
  if (is.null(Mo)) {
    Mo <- N*mean(Mi)
  }
  n <- length(Mi)
  M.bar <- Mo/N
  mu.hat <- sum(taui)/sum(Mi)
  tau.hat <- Mo*mu.hat
  M.bar <- Mo/N
  temp <- (taui - Mi*mu.hat)^2
  varmu.hat <- (1/M.bar^2)*(1-n/N)*sum(temp)/(n*(n-1))
  vartau.hat <- Mo^2*varmu.hat
  Bmu.hat <- 2*sqrt(varmu.hat)
  Btau.hat <- 2*sqrt(vartau.hat)
  Limu <- mu.hat - Bmu.hat; Lsmu <- mu.hat + Bmu.hat
  Litau <- tau.hat - Btau.hat; Lstau <- tau.hat + Btau.hat
  result <- data.frame(Estimation = c(tau.hat, mu.hat),
                       B = c(Btau.hat, Bmu.hat),
                       LI = c(Litau, Limu),
                       LS = c(Lstau, Lsmu))
  rownames(result) <- c("Tau", "Mu")
  S2rcon <- sum(temp)/(n-1)
  return(list(result = result, S2rcon = S2rcon))
}

sample_size <- function(N, sigma2, D) {
  
  #N es la cantidad de conglomerados
  #sigma2 es la estimacion de la varianza
  #D es B^2Mbar^2/4 aunque depende del parametro
  
  n <- N*sigma2/(N*D + sigma2)
  result <- data.frame(nreal = n, ninteger = ceiling(n))
  rownames(result) <- c("ns")
  return(result)
}

#Conglomerados tamaños distintos PPT, propociones y totales

Estimacionesppt <- function(M_i, t_i, n, N, M_0){
  #t_i total por conglomerado
  #N total de conglomerados
  #n número de conglomerados en la muestra
  #M_i Tamaño del iésimo conglomerado
  
  p_i <- M_i/M_0 
  t_pi <- t_i/p_i 
  mu_i <- t_i/M_i 
  t_ppt <- M_0/n*sum(mu_i) 
  mu_ppt <- t_ppt/M_0  
  vart_ppt <- M_0^2/n*sum((mu_i-mu_ppt)^2/(n-1)) 
  varmu_ppt <- 1/n*sum((mu_i-mu_ppt)^2/(n-1)) 
  
  varianzas <- data.frame(T_ppt = vart_ppt, Mu_ppt = varmu_ppt)
  LEEs <- data.frame(T_ppt = 2*sqrt(vart_ppt), Mu_ppt = 2*sqrt(varmu_ppt))
  intervalos <- data.frame(LI = c(t_ppt - 2*sqrt(vart_ppt), mu_ppt - 2*sqrt(varmu_ppt)),
                           LS = c(t_ppt + 2*sqrt(vart_ppt), mu_ppt + 2*sqrt(varmu_ppt)))
  rownames(intervalos) <- c("T_ppt", "Mu_ppt")
  list(Estimaciones = data.frame(T_ppt = t_ppt, Mu_ppt = mu_ppt),
       Mu_i = mu_i, Varianzas = varianzas, LEEs = LEEs, Intervalos = intervalos)
}

samplesizeppt <- function(B, N, n = NULL, mu_i = NULL, mu_ppt = NULL, est = NULL, aprox = T, alpha = NULL){
  #B Límite de error para la estimación
  #N Total de conglomerados 
  #El resto de parámetros se usan si se proporcionan, de lo contrario debe dar el parámetro "est"
  z <- ifelse(aprox, 2, qnorm(1-alpha/2))
  n <- ifelse(is.null(n), (z/B)^2*est, (z/B)^2*sum((mu_i-mu_ppt)^2/(n-1)))
  data.frame(napprox = n, n = ceiling(n))
}

EstimacionesPA <- function(A_i, M_i, N, n, M_0 = NULL){
  #A_i Nro de elementos de interés en el i-ésimo conglomerado
  #M_i Tamaño del i-ésimo conglomerado
  #N Total de conglomerados
  #n número de conglomerados en la muestra 
  if (is.null(M_0)) {
    M_0 <- N*mean(M_i)
  }
  p_con <- sum(A_i)/sum(M_i)
  Mbar <- sum(M_i)/n
  varp_con <- (N-n)/(N*n*Mbar^2)*(S2pcon <- sum((A_i-p_con*M_i)^2)/(n-1))
  A_con <- M_0*p_con
  varA_con <- M_0^2*varp_con
  varianzas <- data.frame(P_con = varp_con, A_con = varA_con)
  LEEs <- data.frame(P_con = 2*sqrt(varp_con), A_con = 2*sqrt(varA_con))
  intervalos <- data.frame(LI = c(p_con - 2*sqrt(varp_con), A_con - 2*sqrt(varA_con)),
                           LS = c(p_con + 2*sqrt(varp_con), A_con + 2*sqrt(varA_con)))
  rownames(intervalos) <- c("P_con", "A_con")
  list(Estimaciones = data.frame(P_con = p_con, A_con = A_con), Varianzas = varianzas,
       LEEs = LEEs, Intervalos = intervalos, S2pcon = S2pcon)
}
