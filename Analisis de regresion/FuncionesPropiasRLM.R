library(tidyverse)
library(MASS)
library(car)
library(lmtest)
library(nortest)
PRESS <- function(linear.model) {
  # calculate the predictive residuals
  pr <- residuals(linear.model) / (1-lm.influence(linear.model)$hat)
  # calculate the PRESS
  PRESS <- sum(pr^2)
  return(PRESS)
}

pred_r_squared <- function(linear.model) {
  # Use anova() to get the sum of squares for the linear model
  lm.anova <- anova(linear.model)
  # Calculate the total sum of squares
  tss <- sum(lm.anova$'Sum Sq')
  # Calculate the predictive R^2
  pred.r.squared <- 1-PRESS(linear.model)/(tss)
  return(pred.r.squared)
}

PruebaFormal.Falta <- function(datos, mod){
  n <- length(datos$x)           
  m <- length(unique(datos$x))   
  datos <- datos %>% 
    group_by(x) %>%
    summarise(f = sum((y - mean(y))^2), .groups = 'drop') %>%
    ungroup()
  SSRES <- anova(mod)$'Sum Sq'[2]
  SSEP <- sum(datos$f)
  SSFA <- SSRES - SSEP
  MSFA <- SSFA/(m - 2)
  MSEP <- SSEP/(n - m)
  F0 <- MSFA/MSEP
  data.frame(SSRES = SSRES, SSEP = SSEP,
            SSFA = SSFA, MSEP = MSEP, MSFA = MSFA,
             F0 = F0, 
             P.Value = pf(q = F0, df1 = m - 2, df2 = n - m, lower.tail = FALSE))
}

Extrapolacion.Oculta <- function(mod = NULL, X = NULL, x0) {
  if(!is.null(mod)){
    hmax <- max(hatvalues(mod))
    if (!("x" %in% names(mod))){
      mod <- update(mod, x = T)
    }
    X <- mod$x
  } else {
    hmax <- max(diag(X %*% solve(t(X) %*% X) %*% t(X)))
  }
  
  if (is.null(dim(x0))){
    x0 <- matrix(x0, nrow = 1)
  }
  h00 <- diag(x0 %*% solve(t(X) %*% X) %*% t(x0))
  data.frame(h00 = h00, hmax = hmax, Extrapolacion = hmax < h00 )
}

Multiple.test <- function(mod, T., r, C){
  if (!("x" %in% names(mod))){
    mod <- update(mod, x = T)
  }
  if (is.null(dim(C))) {
    C <- matrix(C, ncol = 1)
  }
  n <- length(mod$fitted.values)
  p <- mod$rank
  X <- mod$x
  B <- matrix(as.numeric(coef(mod)))
  MSRes <- summary(mod)$sigma^2
  num1 <- (T. %*% B) - C
  num2 <- T. %*% solve(t(X) %*% X) %*% t(T.)
  num <- (t(num1) %*% solve(num2) %*% num1)/r
  den <- MSRes
  F0 <- num/den
  p.value <- pf(q = F0, df1 = r, df2 = n-p, lower.tail = F)
  data.frame(F0 = F0, p_value = p.value, n = n, p = p, SSR.given = num, MSRes = MSRes)
}

SSR <- function(mod){
  mod <- update(mod, x = T, y = T)
  B <- matrix(as.numeric(coef(mod)), ncol = 1)
  X <- mod$x
  Y <- mod$y
  n <- length(Y)
  as.numeric((t(X %*% B) %*% Y) - n*mean(Y)^2)
}

eval.model <- function(mod, values){
  sum(as.numeric(coef(mod))*values)
} 

lambda.optim <- function(mod) {
  bc <- boxcox(mod)
  bc$x[which.max(bc$y)] 
}

anova.montgo <- function(mod) {
   SSR.mod <- SSR(mod)
   SSRes.mod <- anova(mod)$'Sum Sq'[mod$rank]
   SST.mod <- SSR.mod + SSRes.mod
   df <- c(mod$rank - 1, length(fitted(mod)) - mod$rank, length(fitted(mod)) - 1)
   MSQ <- c(SSR.mod/df[1], SSRes.mod/df[2], NA)
   F0 <- c(MSQ[1]/MSQ[2], NA, NA)
   pvalue <- c(pf(q = F0[1], df1 = df[1], df2 = df[2], lower.tail = FALSE), NA,
               NA)
   data.frame(Variacion = c('Regresion', 'Residual', 'Total'),
              Sum.Sq = c(SSR.mod, SSRes.mod, SST.mod), df = df,
              MSQ = MSQ, F0 = F0, p.value = pvalue)
}

#-------------------------------------------------------------------
#EXPERIMENToS PARA REDES NEURONALES Y KAGGLE

# escalamiento <- function(x) {
#   x.escalado <- (x - mean(x))/sqrt(sum((x - mean(x))^2))
# }
# 
# escalamiento2.0 <- function(x) {
#   for (i in 1:(dim(x)[2])) {
#     x[ ,i] <- (x[ ,i] - mean(x[ ,i]))/sqrt(sum((x[,i] - mean(x[ ,i]))^2)) 
#   }
# }
# 
# #Para redes neuronales
# normalize <- function(x) {
#   return((x - min(x)) / (max(x) - min(x)))
# }
# 
# unnormalize <- function(x) {
#   return(x * (max(byke.train$Total) - min(byke.train$Total)) + min(byke.train$Total))
# }
# 
# softplus <- function(x) { log(1 + exp(x)) }
# 
# 
# mapAIC <- function(mod) {
#   mod <- update(mod, x = TRUE)
#   df <- dim(mod$x)[1] - (mod$anova)$'Resid. Df'[length((mod$anova)$'Resid. Df')] + 1
# }








