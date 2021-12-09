################################################################################
#JUANJO
#R tiene implementado el test con la función pero a veces difiere con las diapos
mann_whitney <- function(X, Y, stat = "T1", alternative = "two_sided"){
  #X: Observaciones de la poblacion 1
  #Y: Observaciones de la población 2
  #NOTA: Si se usa el estadístico T, W NO PUEDE SER NULO
  n <- length(X)
  m <- length(Y)
  N <- n + m
  merge_obs <- c(X, Y)
  Ri <- rank(merge_obs)
  Rx <- Ri[1:n]
  Ry <- Ri[(n + 1):N]
  
  Ts <- sum(Rx)
  #Calculos dependen de la selección del estadístico de prueba
  if(stat == "T"){
    Tp <- n*(N + 1) - Ts
    if(alternative == "two_sided"){
      minimum <- min(Ts, Tp)
      num <- minimum - n*(N + 1)/2 + 1/2
      den <- sqrt(n*m*(N + 1)/12)
      proof_stat <- num/den
      pvalue <- 2*pnorm(proof_stat)
    }
    if(alternative == "less"){
      num <- Ts - n*(N + 1)/2 + 1/2
      den <- sqrt(n*m*(N + 1)/12)
      proof_stat <- num/den
      pvalue <- pnorm(proof_stat)
    }
    if(stat == "greater"){
      num <- Tp - n*(N + 1)/2 + 1/2
      den <- sqrt(n*m*(N + 1)/12)
      proof_stat <- num/den
      pvalue <- pnorm(proof_stat)
    }
    return(pvalue)
  }
  #si se va a usar T1
  else{
    num_NC <- Ts - n*(N + 1)/2
    num_C <- Ts - n*(N + 1)/2 - 1/2
    den <- sqrt(n*m/(N*(N - 1)) * sum(Ri^2) - n*m*(N + 1)^2/(4*(N - 1)))
    T1_NC <- num_NC/den
    T1_C <- num_C/den
    if(alternative == "two_sided"){
      pvalue_NC <- 2*min(pnorm(T1_NC), pnorm(T1_NC, lower.tail = F))
      pvalue_C <- 2*min(pnorm(T1_C), pnorm(T1_NC, lower.tail = F))
    }
    if(alternative == "less"){
      pvalue_NC <- pnorm(T1_NC)
      pvalue_C <- pnorm(T1_C)
    }
    if(alternative == "greater"){
      pvalue_NC <- pnorm(T1_NC, lower.tail = F)
      pvalue_C <- pnorm(T1_C, lower.tail = F)
    }
    return(list(Without_Correct_Factor = pvalue_NC,
                With_Correct_Factor = pvalue_C))
  }
}

#INTERVALO DE CONFIANZA PARA DIFERENCIA DE MEDIAS
mean_diffs_confint <- function(X, Y, W){
  #X: Observaciones de la población 1
  #Y: Observaciones de la población 2
  #W: Es el valor crítico de la tabla A7 de Conover o usando aproximación normal
  n <- length(X)
  m <- length(Y)
  r <- n*m
  k <- round(W - n*(n+1)/2)
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
  overall <- list(N = N, Statistic = T2, D_Squared = D2, S_bar = S_bar, 
                  Critical_value = qchisq(alpha, length(populations) - 1, 
                                          lower.tail = F),
                  Ranks_Sum_Squares = Sj,
                  PValue = pvalue,
                  Reject = Reject)
  return(overall)
}

#SIMON
multcomp.var <- function(population, alpha = 0.05) {
  #populations: lista donde cada entrada es un vector con las observaciones
  #de cada población
  #alpha: Nivel de significancia para la región de rechaza (Se fija en 0.05)
  population <- lapply(population, function(x) abs(scale(x, scale = F)))
  nj <- sapply(population, length) #extrayendo los nj
  k <- length(population) #numero de poblaciones
  N <- sum(nj) #calculando n
  temp <- c() #aca se almacenan todos los datos
  for (i in 1:k) {
    temp <- append(temp, population[[i]])
  }
  rankstemp <- rank(temp) #rangos de todos los datos
  ranks <- list()
  start <- 1; end <- nj[1] #primera pareja start, end
  for (i in 1:k) {
    ranks[[i]] <- rankstemp[start:end] #estos pertenecen a la poblacion i
    if (i < k) {
      start <- end +1; end <- end + nj[i+1] #se actualizan al final, solo si i < k
    }
  }
  Sjtemp <- lapply(ranks, function(x) x^2) #se elevan los rangos al cuadrado por muestra
  Sbar <- sum(Sj <- sapply(Sjtemp, sum))/N #se calcula s barra
  D.sq <- (sum(rankstemp^4) - N*Sbar^2)/(N-1) #se calcula d barra cuadrado
  T.2 <- (sum(Sj^2/nj) - N*Sbar^2)/D.sq #se calcula el estadistico de prueba
  result <- list(N = N, Sbar = Sbar, Dsq = D.sq, T2 = T.2,
                 p.value = pchisq(T.2, k-1, lower.tail = F),
                 reject = pchisq(T.2, k-1, lower.tail = F) < alpha,
                 critical.value = qchisq(alpha, k - 1, lower.tail = F))
  if (result$reject) {
    idx <- combn(1:k, 2)
    reps <- ncol(idx)
    aux <- rep(0, reps)
    comps <- data.frame(Population = aux, diffs = aux, direction = as.character(aux), critic.value = aux)
    for (i in 1:reps) {
      idxtemp <- idx[,i]
      first <- idxtemp[1]
      second <- idxtemp[2]
      comps[i, 1] <- paste(first, second, sep = "-")
      comps[i, 2] <- abs(Sj[first]/nj[first] - Sj[second]/nj[second])
      comps[i, 4] <- qt(alpha/2, N - ncol(idx), lower.tail = F)*sqrt(D.sq*(N-1-T.2)/(N-k))*sqrt(1/nj[first]+1/nj[second])
      comps[i, 3] <- ifelse(comps[i,2] > comps[i,4], ">", "<")
    }
    return(list(result, comps = comps))
  }
  return(result)
}

#Numero de concordantes y discordantes
conc.dis <- function(x, y) {
  library(tidyverse)
  dfaux <- data.frame(x, y)
  dfaux <- dfaux %>%
    arrange(x, y)
  x <- dfaux$x
  y <- dfaux$y
  rm(dfaux)
  n <- length(x)
  df <- data.frame(pair = rep("", n-1), Conc = rep(0, n-1), Disc = rep(0, n-1))
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      if (x[i] != x[j]) {
        if((y[i] - y[j])/(x[i] - x[j]) > 0) {
          df[i, 2] <- df[i, 2] + 1
        } else if ((y[i] - y[j])/(x[i] - x[j]) < 0) {
          df[i, 3] <- df[i, 3] + 1
        } else {
          df[i, 2] <- df[i, 2] + 0.5
          df[i, 3] <- df[i, 3] + 0.5
        }
      } else {
        next
      }
    }
    df[i, 1] <- paste("(", x[i], ", ", y[i], ")", sep = "")
  }
  return(list(df = df, Nc = sum(df[,2]), Nd = sum(df[,3])))
}
