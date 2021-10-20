binomial.exact.test <- function(x, n, p, kind = "two_sided", alpha = NULL) {
  #x es el numero de exitos
  #n es el tamano de muestra
  #p es la probabilidad a contrastar
  #kind es left, right o two_sided (default)
  #alpha es el nivel de significancia, si es de dos colas debe ser un vector
  if (!is.null(alpha)) {
    if (kind == "left") {
      for (i in 0:n) {
        if (pbinom(i, n, p) > alpha) {
          return(list(lower.limit = i-1, alpha.real = pbinom(i-1, n, p)))
        }
      }
    } else if (kind == "right") {
      for (i in n:0) {
        if (pbinom(i, n, p, lower.tail = F) > alpha) {
          return(list(upper.limit = i+1, alpha.real = pbinom(i+1, n, p, lower.tail = F)))
        }
      }
    } else {
      lower <- 0; upper <- 0
      for (i in 0:n) {
        if (pbinom(i, n, p) > alpha[1]) {
          lower <-  i-1
          break
        }
      }
        for (i in n:0) {
          if (pbinom(i, n, p, lower.tail = F) > alpha[2]) {
            upper <- i+1 
            break
          }
        }
      return(list(lower.limit = lower, upper.limit = upper, alpha.real = pbinom(lower, n, p) + pbinom(upper, n, p,lower.tail = F)))
    }
    
  } else {
    if (kind == "left") {
      p.value <- pbinom(x, n, p)
    } else if (kind == "right") {
      p.value <- pbinom(x-1, n, p, lower.tail = F)
    } else {
      p.value <- 2*min(pbinom(x, n, p), pbinom(x-1, n, p, lower.tail = F))
    }
    return(list(p.value = p.value))
  }
}

binomial.asintotic.test <- function(x, n, p, kind = "two_sided", alpha= NULL) {
  #funciona igual que el de arriba, con n > 20
  #alpha en este caso es un solo valor siempre sin importar la naturaleza de la prueba
  if (is.null(alpha)) {
    den <- sqrt(n*p*(1-p))
    numpar <- x - n*p
    if (kind == "left") {
      p.value <- pnorm((numpar + 0.5)/den)
    } else if (kind == "right") {
      p.value <- pnorm((numpar - 0.5)/den, lower.tail = F)
    } else {
      p.value <- 2*min(pnorm((numpar + 0.5)/den), pnorm((numpar - 0.5)/den, lower.tail = F))
    }
    return(list(p.value = p.value))
  } else {
    fact <- sqrt(n*p*(1-p))
    if (kind == "left") {
      return(list(lower.limit = n*p + qnorm(alpha)*fact))
    } else if (kind == "right") {
      return(list(upper.limit = n*p + qnorm(alpha, lower.tail = F)*fact))
    } else {
      return(list(lower.limit = n*p + qnorm(alpha/2)*fact,
                  upper.limit = n*p + qnorm(alpha/2, lower.tail = F)*fact))
    }
  }
}


#TEST CUANTIL
exact_quantile_test <- function(data, p, q, kind = "two_sided", alpha = NULL){
  #data: Vector de datos
  #p: probabilidad p*
  #q: cuantil x*
  #kind: Tipo de prueba (two_sided, left, right)
  #alpha: nivel de significancia de la prueba; si no es dado se calcula el exacto.
  #Si es prueba de dos colas, alpha es un vector de dos entradas siendo la primera
  #el area a izquierda y la segunda el area a derecha.
  #Nota: Para el test de dos colas no se muestra información para calcular el 
  #valor-p con la distribucion exacta por lo que cuando alpha es nulo solo se
  #calculan los valores-p en pruebas de cola izquierda o derecha
  n <- length(data)
  T1 <- sum(data <= q); T2 <- sum(data < q) #necesario para la RR
  if(!is.null(alpha)){
    if(kind == "two_sided"){
      for(i in 0:n){
        if(pbinom(i, n, p) > alpha[1]){
          t1 <- i - 1
          alpha1 <- pbinom(t1, n, p)
          break
        }
      }
      for(i in n:0){
        if(pbinom(i, n, p, lower.tail = F) > alpha[2]){
          t2 <- i + 1
          alpha2 <- pbinom(t2, n, p, lower.tail = F)
          break
        }
      }
      actual_alpha <- alpha1 + alpha2
      if(T1 <= t1 | T2 > t2){
        Reject <- T
      }
      else{
        Reject <- F
      }
      info <- data.frame(T_Statistics = c(T1, T2),
                         t_Quantiles = c(t1, t2),
                         alphas = c(alpha1, alpha2))
      overall <- list(Info = info,
                      actual_alpha = actual_alpha,
                      Reject = Reject)
      return(overall)
    }
    if(kind == "left") {
      for(i in 0:n) {
        if(pbinom(i, n, p) > alpha){
          t1 <- i - 1
          actual_alpha <- pbinom(t1, n, p)
          break
        }
      }
      if(T1 <= t1){
        Reject <- T
      }
      else{
        Reject <- F
      }
      info <- data.frame(T_Statistic = T1,
                         t_Quantile = t1)
      overall <- list(Info = info,
                      actual_alpha <- actual_alpha,
                      Reject = Reject)
      return(overall)
    }
    if(kind == "right"){
      for(i in n:0){
        if(pbinom(i, n, p, lower.tail = F)){
          t2 <- i + 1
          actual_alpha <- pbinom(t2, n, p, lower.tail = F)
          break
        }
      }
      if(T2 > t2){
        Reject <- T
      }
      else{
        Reject <- F
      }
      info <- data.frame(T_Statistic = T2, 
                         t_Quantile = t2)
      overall <- list(Info = info,
                      actual_alpha = actual_alpha,
                      Reject = Reject)
      return(overall)
    }
  }
  else{
    if(kind == "left"){
      pvalue <- pbinom(T1, n, p)
      overall <- list(pvalue = pvalue)
      return(overall)
    }
    else if(kind == "right"){
      pvalue <- pbinom(T2 - 1, n, p, lower.tail = F)
      overall <- list(pvalue = pvalue)
      return(overall)
    }
    else{
      return(print("Not suported P-Value for a two sided test"))
    }
  }
}

asintotic_quantile_test <- function(p, q, data = NULL, n = NULL, T1 = NULL, T2 = NULL,
                                    kind = "two_sided", alpha = NULL){
  #data: Vector de datos
  #Si data es nulo, se debe pasar el estadistico de prueba y el tamaño de muestra
  #p: probabilidad p*
  #q: cuantil x*
  #kind: Tipo de prueba de las colas (two-sided, left, right)
  #alpha: Nivel de significancia de la prueba, si es nulo se usa region de rechazo
  if(!is.null(data)){
    n <- length(data)
    T1 <- sum(data <= q); T2 <- sum(data < q)
    if(!is.null(alpha)){
      if(kind == "two_sided"){
        t1 <- n*p + qnorm(alpha/2) * sqrt(n*p*(1 - p))
        t2 <- n*p - qnorm(alpha/2, lower.tail = F) * sqrt(n*p*(1 - p))
        if(T1 <= t1 | T2 > t2){
          Reject <- T
        }
        else{
          Reject <- F
        }
        info <- data.frame(T_Statistics = c(T1, T2),
                           t_Quantiles = c(t1, t2))
        overall <- list(Info = info, 
                        Reject = Reject)
        return(overall)
      }
      if(kind == "left"){
        t1 <- n*p + qnorm(alpha) * sqrt(n*p*(1 - p))
        if(T1 <= t1){
          Reject <- T
        }
        else{
          Reject <- F
        }
        info <- data.frame(T1_Statistic = T1,
                           t1_Quantile = t1)
        overall <- list(Info = info, 
                        Reject = Reject)
        return(overall)
      }
      if(kind == "right"){
        t2 <- n*p - qnorm(alpha) * sqrt(n*p*(1 - p))
        if(T2 > t2){
          Reject <- T
        }
        else{
          Reject <- F
        }
      }
    }
    else{
      den <- sqrt(n*p(1 -p))
      proof_stat1 <- (T1 - n*p + 0.5)/den
      proof_stat2 <- (T2 - n*p - 0.5)/den
      if(kind == "two_sided"){
        P1 <- pnorm(proof_stat1)
        P2 <- pnorm(proof_stat2, lower.tail = F)
        pvalue <- 2 * min(P1, P2)
      }
      else if(kind == "left"){
        pvalue <- pnorm(proof_stat1)
      }
      else{
        pvalue <- pnorm(proof_stat2, lower.tail = F)
      }
      return(list(pvalue = pvalue))
    }
  }###############################
  else{
    if(!is.null(alpha)){
      if(kind == "two_sided"){
        t1 <- n*p + qnorm(alpha/2) * sqrt(n*p*(1 - p))
        t2 <- n*p - qnorm(alpha/2, lower.tail = F) * sqrt(n*p*(1 - p))
        if(T1 <= t1 | T2 > t2){
          Reject <- T
        }
        else{
          Reject <- F
        }
        info <- data.frame(T_Statistics = c(T1, T2),
                           t_Quantiles = c(t1, t2))
        overall <- list(Info = info, 
                        Reject = Reject)
        return(overall)
      }
      if(kind == "left"){
        t1 <- n*p + qnorm(alpha) * sqrt(n*p*(1 - p))
        if(T1 <= t1){
          Reject <- T
        }
        else{
          Reject <- F
        }
        info <- data.frame(T1_Statistic = T1,
                           t1_Quantile = t1)
        overall <- list(Info = info, 
                        Reject = Reject)
        return(overall)
      }
      if(kind == "right"){
        t2 <- n*p - qnorm(alpha) * sqrt(n*p*(1 - p))
        if(T2 > t2){
          Reject <- T
        }
        else{
          Reject <- F
        }
      }
    }
    else{
      den <- sqrt(n*p*(1 -p))
      proof_stat1 <- (T1 - n*p + 0.5)/den
      proof_stat2 <- (T2 - n*p - 0.5)/den
      if(kind == "two_sided"){
        P1 <- pnorm(proof_stat1)
        P2 <- pnorm(proof_stat2, lower.tail = F)
        pvalue <- 2 * min(P1, P2)
      }
      else if(kind == "left"){
        pvalue <- pnorm(proof_stat1)
      }
      else{
        pvalue <- pnorm(proof_stat2, lower.tail = F)
      }
      return(list(pvalue = pvalue))
    }
  }
}

#TEST DEL SIGNO
sign_test <- function(alpha, n = NULL, data = NULL, T_stat = NULL, kind = "two_sided"){
  #data: bd con las observaciones de las dos v.a (suponiendo Xi es la columna 1
  # y Yi es la columna 2)
  #Notas: - Como p = 0.5, la distribucion es simetrica
  #- De la data se puede deducir el n y el estadistico de prueba, si data es nulo
  #se deben especificar dichos argumentos
  #- No se da información sobre valores-p usando la distribucion exacta por lo
  #que se asume que la muestra es lo suficientemente grande cuando se usan 
  #valores-p (alpha nulo)
  if(!is.null(data)){
    n <- dim(data)[1]
    T_stat <- sum(data[, 1] < data[, 2])
    if(n <= 20){
      #En este caso se usa la distribucion exacta
      if(kind == "two_sided"){
        for(i in 0:n){
          if(pbinom(i, n, 0.5) > alpha/2){
            t <- i - 1
            break
          }
        }
        if(T_stat <= t | T_stat >= n - t){
          Reject <- T
        }
        else{
          Reject <- F
        }
        pvalue <- 2 * min(pbinom(T_stat, n, 0.5, lower.tail = F),
                          pbinom(T_stat, n, 0.5))
        actual_alpha <- 2 * pbinom(t, n, 0.5)
      }
      else if(kind == "right"){
        for(i in 0:n){
          if(pbinom(i, n, 0.5) > alpha){
            t <- i - 1
            break
          }
        }
        if(T_stat >= n - t){
          Reject <- T
        }
        else{
          Reject <- F
        }
        pvalue <- pbinom(T_stat, n, 0.5, lower.tail = F)
        actual_alpha <- pbinom(t, n, 0.5)
      }
      else {
        for(i in 0:n){
          if(pbinom(i, n, 0.5) > alpha){
            t <- i - 1
            break
          }
        }
        if(T_stat <= t){
          Reject <- T
        }
        else{
          Reject <- F
        }
        pvalue <- pbinom(T_stat, n, 0.5)
        actual_alpha <- pbinom(t, n, 0.5)
      }
      overall <- list(T_statistic = T_stat,
                      t_quantile = t,
                      Reject = Reject,
                      pvalue = pvalue)
      return(overall)
    }
    #n > 20
    else{
      t <- 0.5 * (n + qnorm(alpha/2) * sqrt(n))
      if(kind == "two_sided"){
        num1 <- T_stat - n*p + 0.5
        num2 <- T_stat - n*p - 0.5
        den <- sqrt(n*p*(1 - p))
        pvalue <- 2 * min(pnorm(num1/den),
                          pnorm(num2/den, lower.tail = F))
      }
      else if(kind == "right"){
        num2 <- T_stat - n*p - 0.5
        den <- sqrt(n*p*(1 - p))
        pvalue <- pnorm(num2/den, lower.tail = F)
      }
      else {
        num1 <- T_stat - n*p + 0.5
        den <- sqrt(n*p*(1 - p))
        pvalue <- pnorm(num2/den)
      }
      overall <- list(T_statistic = T_stat,
                      t_quantile = t,
                      pvalue = pvalue)
      return(overall)
    }
  }
  #data nula
  else{
    if(n <= 20){
      #En este caso se usa la distribucion exacta
      if(kind == "two_sided"){
        for(i in 0:n){
          if(pbinom(i, n, 0.5) > alpha/2){
            t <- i - 1
            break
          }
        }
        if(T_stat <= t | T_stat >= n - t){
          Reject <- T
        }
        else{
          Reject <- F
        }
        pvalue <- 2 * min(pbinom(T_stat, n, 0.5, lower.tail = F),
                          pbinom(T_stat, n, 0.5))
        actual_alpha <- 2 * pbinom(t, n, 0.5)
      }
      else if(kind == "right"){
        for(i in 0:n){
          if(pbinom(i, n, 0.5) > alpha){
            t <- i - 1
            break
          }
        }
        if(T_stat >= n - t){
          Reject <- T
        }
        else{
          Reject <- F
        }
        pvalue <- pbinom(T_stat, n, 0.5, lower.tail = F)
        actual_alpha <- pbinom(t, n, 0.5)
      }
      else {
        for(i in 0:n){
          if(pbinom(i, n, 0.5) > alpha){
            t <- i - 1
            break
          }
        }
        if(T_stat <= t){
          Reject <- T
        }
        else{
          Reject <- F
        }
        pvalue <- pbinom(T_stat, n, 0.5)
        actual_alpha <- pbinom(t, n, 0.5)
      }
      overall <- list(T_statistic = T_stat,
                      t_quantile = t,
                      Reject = Reject,
                      pvalue = pvalue)
      return(overall)
    }
    #n > 20
    else{
      t <- 0.5 * (n + qnorm(alpha/2) * sqrt(n))
      if(kind == "two_sided"){
        num1 <- T_stat - n*p + 0.5
        num2 <- T_stat - n*p - 0.5
        den <- sqrt(n*p*(1 - p))
        pvalue <- 2 * min(pnorm(num1/den),
                          pnorm(num2/den, lower.tail = F))
      }
      else if(kind == "right"){
        num2 <- T_stat - n*p - 0.5
        den <- sqrt(n*p*(1 - p))
        pvalue <- pnorm(num2/den, lower.tail = F)
      }
      else {
        num1 <- T_stat - n*p + 0.5
        den <- sqrt(n*p*(1 - p))
        pvalue <- pnorm(num2/den)
      }
      overall <- list(T_statistic = T_stat,
                      t_quantile = t,
                      pvalue = pvalue)
      return(overall)
    }
  }
}


