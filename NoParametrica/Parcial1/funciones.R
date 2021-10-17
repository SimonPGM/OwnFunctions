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

asintotic_quantile_test <- function(data, p, q, kind = "two_sided", alpha = NULL){
  #data: Vector de datos
  #p: probabilidad p*
  #q: cuantil x*
  #kind: Tipo de prueba de las colas (two-sided, left, right)
  #alpha: Nivel de significancia de la prueba, si es nulo se usa region de rechazo
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
}
