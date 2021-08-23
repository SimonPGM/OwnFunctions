#Conglomerados con tamaños iguales
mu_tau_con <- function(tau_i, M, N, n, approx = T, alpha = NULL){
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
  
  #Errores estandars
  if(approx){
    B_mu_c <- 2 * sqrt(var_mu_c_hat)
    B_mu_hat <- 2 * sqrt(var_mu_hat)
    B_tau_hat <- 2 * sqrt(var_tau_hat)
  }
  else{
    B_mu_c <- qnorm(alpha/2, lower.tail = F) * sqrt(var_mu_c_hat)
    B_mu_hat <- qnorm(alpha/2, lower.tail = F) * sqrt(var_mu_hat)
    B_tau_hat <- qnorm(alpha/2, lower.tail = F) * sqrt(var_tau_hat)
  }
  
  #Limites de error de estimacion
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
                  LEE = LEE,
                  ICs = IC)
  return(overall)
}

#Conglomerados tamaños distintos MAS y Razon


#Conglomerados tamaños distintos PPT, propociones y totales