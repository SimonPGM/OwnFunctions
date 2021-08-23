samples <- function(sdata, n){
  sdata <- sort(sdata)
  distdatos <- prop.table(table(sdata)) #distribucion de los datos
  mu <- sum(distdatos*unique(sdata)) #media poblacional
  sigma2 <- sum(distdatos*(unique(sdata)-mu)^2) #varianza poblacional
  combinations <- combn(sdata, n) #todas las muestras posibles
  ncombinations <- length(combinations) #numero de combinaciones posibles
  ybarvalues <- sort(apply(combinations, 2, mean)) #valores posibles de la media muestral
  distybar <- prop.table(table(ybarvalues)) #distribucion de la media muestral
  sigma2ybar <- sum(distybar*(unique(ybarvalues)-mu)^2) #varianza de la media muestral
  s2values <-  sort(apply(combinations, 2, var)) #valores de s2
  s2dist <- prop.table(table(s2values)) #disitribucion varianza muestral
  s2esp <- sum(s2dist*unique(s2values)) #valor esperado varianza muestral
  list(datos = sdata, distdatos = distdatos, mu_pob = mu, var_pob = sigma2,
       combinaciones = combinations, n_combinaciones = ncombinations, valores_ybarra = ybarvalues,
       distybar = distybar, varybar = sigma2ybar, valoress2 = s2values, dists2 = s2dist, expecteds2 = s2esp)
}

confint_p_A <- function(p, n, N, alpha) {
  status = T
  if (p > 0.5){
    if (n <= 30){
      status <- F
    }
  }
  else if (p < 0.005){
    if (n <= 1400){
      status <- F
    }
  }
  else{
    ps <- c(0.5, 0.4, 0.3, 0.2, 0.1, 0.05)
    ns <- c(30, 50, 80, 200, 600, 1400)
    difi <- which.min(abs(p-ps))
    if (n <= ns[difi]){
      status <- F
    }
  }
  b <- ifelse(status, qnorm(1-alpha/2), qt(1-alpha/2, n-1))*sqrt((1-n/N)*p*(1-p)/(n-1))
  list(liminfp = p - b, limsupp = p + b, liminfA = N*(p - b), limsupA = N*(p + b), A = N*p, norm = status)
}

confint_mu_tau <- function(datos = NULL, mu = NULL, s2 = NULL, n, N, alpha) {
  status <- n >30
  mu <- ifelse(is.null(datos), mu,mean(datos))
  s2 <- ifelse(is.null(datos), s2, var(datos))
  b <- ifelse( status, qnorm(1-alpha/2), qt(1-alpha/2, n-1))*sqrt((1-n/N)*s2/n)
  list(mu = mu, tau = N*mu, liminfmu = mu - b, limsupmu = mu + b, liminftau = N*(mu - b),
       limsuptau = N*(mu + b), norm = status)
}

sample_size <- function(N, B, sigma2, alpha, total = F){
  n0 <- ifelse(total, N^2,1)*qnorm(1-alpha/2)^2*sigma2/B^2
  n <- ceiling(1/(1/N + (1-1/N)/n0))
  n
}