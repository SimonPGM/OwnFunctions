################################################################################
#JUANJO

#INTERVALO DE CONFIANZA PARA DIFERENCIA DE MEDIAS
mean_diffs_confint <- function(X, Y, W){
  #X: Observaciones de la población 1
  #Y: Observaciones de la población 2
  #W: Es el valor crítico de la tabla A7 de Conover o usando aproximación normal
  n <- length(X)
  m <- length(Y)
  r <- n*m
  k <- W - n*(n+1)/2
  aux <- rep(X, each = length(Y))
  all_diffs <- aux - Y
  all_diffs_sort <- sort(all_diffs)
  IC <- all_diffs_sort[c(k, r - k + 1)]
  IC
}

#TEST PARA COMPARAR VARIANZAS DE VARIAS POBLACIONES
comp_multiple_vars <- function(populations, alpha = 0.05){
  #populations: lista donde cada entrada es un vector con las observaciones
  #de cada población
  #alpha: Nivel de significancia para la región de rechaza (Se fija en 0.05)
  
  #Tamaño de la j-esima muestra
  nj <- sapply(populations, length)
  N <- sum(nj)
  ui <- list()
  for(i in 1:length(populations)){
    ui[[i]] <- abs(populations[[i]] - mean(populations[[i]]))
  }
  
  #Agrupando todo en un solo vector
  merge_obs <- c()
  for(i in 1:length(ui)){
    merge_obs <- append(merge_obs, ui[[i]])
  }
  #Rangos de todas las observaciones
  ranks <- rank(merge_obs)
  
  #Obteniendo los rangos de cada muestra
  Ri <- list()
  from <- 1
  to <- length(ui[[1]])
  for(i in 1:length(ui)){
    Ri[[i]] <- ranks[from:to]
    from <- from + length(ui[[i]])
    if(i < length(ui)){
      to <- to + length(ui[[i + 1]]) 
    }
  }
  #Suma de cuadrados de los rangos
  Sj <- c()
  for(i in 1:length(Ri)){
    Sj[i] <- sum(Ri[[i]]^2)
  }
  S_bar <- (1/N)*sum(Sj)
  
  #Suma de rangos a la cuarta potencia
  Ri4 <- c()
  for(i in 1:length(Ri)){
    Ri4[i] <- sum(Ri[[i]]^4)
  }
  D2 <- (1/(N - 1))*(sum(Ri4) - N*S_bar^2)
  
  #Estadístico de prueba
  T2 <- (1/D2)*(sum((Sj^2)/nj) - N*S_bar^2)
  
  #Región de rechazo
  Reject <- T2 > qchisq(alpha, length(populations) - 1, lower.tail = F)
  
  #Valor-p
  pvalue <- pchisq(T2, length(populations) - 1, lower.tail = F)
  
  #Resumen
  overall <- list(Statistic = T2,
                  Ranks_Sum_Squares = Sj,
                  PValue = pvalue,
                  Reject = Reject)
  return(overall)
}
