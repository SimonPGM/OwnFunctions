#Conglomerados con tamaños iguales


#Conglomerados tamaños distintos MAS y Razon


#Conglomerados tamaños distintos PPT, propociones y totales

Estimacionesppt <- function(M_i, t_i, n, N, M_0 = NULL){
  #t_i total por conglomerado
  #N total de conglomerados
  #n número de conglomerados en la muestra
  #M_i Tamaño del iésimo conglomerado
  if (is.null(M_0)){
    M_0 <- sum(M_i)
  }
  p_i <- M_i/M_0 
  t_pi <- t_i/p_i 
  mu_i <- t_i/M_i 
  t_ppt <- M_0/n*sum(mu_i) 
  mu_ppt <- t_ppt/M_0  
  vart_ppt <- M_0^2/n*sum((mu_i-mu_ppt)^2/(n-1)) 
  varmu_ppt <- 1/n*sum((mu_i-mu_ppt)^2/(n-1)) 
  intervalos <- data.frame(LI = c(t_ppt - 2*sqrt(vart_ppt), mu_ppt - 2*sqrt(varmu_ppt)),
                           LS = c(t_ppt + 2*sqrt(vart_ppt), mu_ppt + 2*sqrt(varmu_ppt)))
  rownames(intervalos) <- c("T_ppt", "Mu_ppt")
  list(Estimaciones = data.frame(T_ppt = t_ppt, Mu_ppt = mu_ppt),
       Mu_i = mu_i, Intervalos = intervalos)
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
  if (is.null(M_0)){
    M_0 <- sum(M_i)
  }
  p_con <- sum(Ai)/sum(Mi)
  Mbar <- 1/N*sum(M_i)
  varp_con <- (N-n)/(N*n*Mbar^2)*sum((A_i-p_con*sum(M_i))^2)/(n-1)
  A_con <- M_0*p_con
  varA_con <- M_0^2*varp_con
  intervalos <- data.frame(LI = c(p_con - 2*sqrt(varp_con), A_con - 2*sqrt(varA_con)),
                           LS = c(p_cont + 2*sqrt(varp_con), A_con + 2*sqrt(varA_con)))
  rownames(intervalos) <- c("P_con", "A_con")
  ist(Estimaciones = data.frame(P_con = p_con, A_con = A_con), Intervalos = intervalos)
}
