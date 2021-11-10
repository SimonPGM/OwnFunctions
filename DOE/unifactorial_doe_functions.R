#Funcion para calcular la tabla anova de efectos aleatorios
anova.random.effects <- function(modelo){
  #modelo : A aov class object
  partial.anova <- summary(modelo)
  names(partial.anova) <- c("factor", "error")
  SSA <- round(as.data.frame.list(partial.anova$factor)[1,2],4)
  SSE <- round(as.data.frame.list(partial.anova$error)[1,2],4)
  MSA <- round(as.data.frame.list(partial.anova$factor)[1,3],4)
  MSE <- round(as.data.frame.list(partial.anova$error)[1,3],4)
  F_0 <- round(MSA/MSE,4)
  dffac <- as.numeric(as.data.frame.list(partial.anova$factor)[1])
  dferr <- as.numeric(as.data.frame.list(partial.anova$error)[1])
  P <- round(pf(F_0, dffac, dferr, lower.tail = F),4)
  
  anova.table <- data.frame(Fuente = c(names(modelo)[2], "Error"), 
                            Df = c(dffac, dferr),
                            "Sum sq" = c(SSA, SSE),
                            "Mean sq" = c(MSA, MSE),
                            "F value" = c(F_0, " "),
                            "Pr(>F)" = c(P, " "))
  return(anova.table)
}

#Hace estimaciones puntuales e intervalos de confianza para los parametros de 
#interés en un DCA unifactorial de efectos aleatorios
estimate_random_efects_anova <- function(anova_table, ni, a, response, gamma = 0.05){
  #anova_table: se calcula con la funcion anova.random.effects
  # a: numero de tratamientos
  #ni: observaciones del tratamiento i
  #response: vector con los valores de la variable respuesta
  SSA <- anova_table$Sum.sq[1]
  SSE <- anova_table$Sum.sq[2]
  MSA <- anova_table$Mean.sq[1]
  MSE <- anova_table$Mean.sq[2] #Estimacion de sigma^2
  if(length(unique(ni) == 1)){
    c = ni[1]
    N <- a * c
  }
  else{
    N <- sum(ni)
    c = (1/(N*(a - 1))) * (N^2 - sum(ni^2))
  }
  #Estimaciones puntuales
  sigma2_alpha <- (MSA - MSE)/c
  mu_hat <- (1/N) * sum(response) 
  
  #Razon de componentes de varianzas sigma2_a/sigma2
  vars_rate <- (MSA - MSE)/(c * MSE) 
  #Proporción de var total debida al factor sigma2_a/(sigma2a + sigma2)
  prop_var_total <- (MSA - MSE)/((c - 1) * MSE + MSA)
  
  #Intervalo para sigma^2
  LI_sigma2 <- SSE/qchisq(gamma/2, N - a, lower.tail = F)
  LS_sigma2 <- SSE/qchisq(gamma/2, N - a)
  IC_sigma2 <- c(LI_sigma2, LS_sigma2)
  
  #Intervalo para sigma^2_alpha
  num_sigma2_alpha <- (MSA - MSE)^2
  den_sigma2_alpha <- (MSA^2)/(a - 1) + (MSE^2)/(N - a)
  kappa <- num_sigma2_alpha/den_sigma2_alpha
  LI_sigma2_alpha <- kappa*sigma2_alpha/qchisq(gamma/2, kappa, lower.tail = F)
  LS_sigma2_alpha <- kappa*sigma2_alpha/qchisq(gamma/2, kappa)
  IC_sigma2_alpha <- c(LI_sigma2_alpha, LS_sigma2_alpha)
  
  #Intervalo para mu_hat
  a1 <- sum(ni^2/(c*N^2))
  a2 <- (1/N - sum(ni^2/(c*N^2))) 
  S2y <- a1 * MSA + a2 * MSE
  num_mu <- (a1*MSA + a2*MSE)^2
  den_mu <- ((a1*MSA)^2)/(a - 1) + ((a2*MSE)^2)/(N - a)
  nu_star <- num_mu/den_mu
  LI_mu <- mu_hat - qt(gamma/2, nu_star, lower.tail = F) * sqrt(S2y)
  LS_mu <- mu_hat + qt(gamma/2, nu_star, lower.tail = F) * sqrt(S2y)
  IC_mu <- c(LI_mu, LS_mu)
  
  #Intervalo para sigma2_a/sigma2
  Fo <- as.numeric(anova_table$F.value[1])
  L <- (1/c)*((Fo/qf(gamma/2, a - 1, N - a, lower.tail = F)) - 1)
  U <- (1/c)*((Fo/qf(gamma/2, a - 1, N - a)) - 1)
  LI_vars_rate <- L
  LS_vars_rate <- U
  IC_vars_rate <- c(LI_vars_rate, LS_vars_rate)
  
  #Intervalo para sigma2_a/(sigma2_a + sigma2)
  LI_prop_var_total <- (L/(1 + L))
  LS_prop_var_total <- (U/(1 + U))
  IC_prop_var_total <- c(LI_prop_var_total, LS_prop_var_total)
  
  #Resumen
  estimations <- data.frame(Sigma_Squared = MSE, 
                            Sigma_Subalpha_Squared = sigma2_alpha,
                            Mu = mu_hat, 
                            Vars_Rate = vars_rate,
                            Prop_Var_Total_Factor = prop_var_total)
  
  ICs <- data.frame(Sigma_Squared = IC_sigma2,
                    Sigma_Subalpha_Squared = IC_sigma2_alpha,
                    Mu = IC_mu,
                    Var_Rate = IC_vars_rate,
                    Prop_Var_Total_Factor = IC_prop_var_total)
  rownames(ICs) <- c("Lower", "Upper")
  overall <- list(Estimations = estimations, 
                  Confidence_Intervals = t(ICs))
  return(overall)
}

#Chequea supuestos para el análisis de residuales
random.effects.residuals <- function(modelo, datos){
  #modelo: lmer object
  #datos: data frame with 2 columns, response and factor respectly
  plot(fitted.values(modelo), scale(residuals(modelo)), 
       main = "Chequeo Homoscedasticidad")
  
  model.residuals <- scale(residuals(modelo))
  
  qqnorm(model.residuals, main = "Chequeo Normalidad de los errores")
  qqline(model.residuals)
  
  sh.residuals <- shapiro.test(model.residuals)
  
  random.effects <- scale(sapply(split(datos, datos[,2]), function(X) mean(X[,1])))
  
  qqnorm(random.effects, main = "Chequeo Normalidad de los efectos aleatorios")
  qqline(random.effects)
  
  sh.effects <- shapiro.test(random.effects)
  
  data.frame(Data = c(sh.residuals$data.name,sh.effects$data.name),
             Pvalue = c(sh.residuals$p.value, sh.effects$p.value))
}

