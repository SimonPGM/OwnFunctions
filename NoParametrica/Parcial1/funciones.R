binomial.exact.test <- function(x, n, p, kind, alpha = NULL) {
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
          return(list(upper.limit = i-1, alpha.real = pbinom(i-1, n, p, lower.tail = F)))
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
            upper <- i-1
            break
          }
        }
      return(list(lower.limit = lower, upper.limit = upper, alpha.real = pbinom(lower, n, p) + pbinom(upper, n, p)))
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

binomial.asintotic.test <- function(x, n, p, kind, alpha= NULL) {
  #funciona igual que el de arriba, con n > 20
  #alpha en este caso es un solo valor siempre sin importar la naturaleza de la prueba
  if (is.null(alpha)) {
    den <- sqrt(n*p*(1-p))
    numpar <- x - n*p
    if (kind == "left") {
      p.value <- pnorm((numpar + 0.5)/den)
    } else if (kind == "right") {
      p.value <- pnorm((numpar - 0.5)/den)
    } else {
      p.value <- 2*min(pnorm((numpar + 0.5)/den), pnorm((numpar - 0.5)/den))
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