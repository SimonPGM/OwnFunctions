#Conglomerados con tamaños iguales


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
  colnames(result) <- c("Tau", "Mu_c", "Mu")
  return(result)
}

tau_muc_mu_r <- function(Mi, taui, N, Mo) {
  
  #Mi tamaño de los conglomerados
  #Taui totales de los conglomerados
  #N total de conglomerados
  #n número de coglomerados seleccionados en la muestra
  #Mo es el total de unidades muestrales en la poblacion
  
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
  colnames(result) <- c("Tau", "Mu")
  return(result)
}

sample_size <- function(N, sigma2, D) {
  
  #N es la cantidad de conglomerados
  #sigma2 es la estimacion de la varianza
  #D es B^2Mbar^2/4
  
  n <- N*sigma2/(N*D + sigma2)
  result <- data.frame(nreal = n, ninteger = ceiling(n))
  colnames(result) <- c("ns")
  return(result)
}

#Conglomerados tamaños distintos PPT, propociones y totales