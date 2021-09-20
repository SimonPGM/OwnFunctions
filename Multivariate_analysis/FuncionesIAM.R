partitions <- function(mu, sigma, p1i, p2i, row_counter = 1, column_counter = 1){
  #mu is the means vector
  #sigma is the variance-covariance matrix
  #pki is the vector of indexes of the k-th partition with k = 1, 2
  
  mu1 <- matrix(mu[p1i], ncol = 1, byrow = T)
  mu2 <- matrix(mu[p2i], ncol = 1, byrow = T)
  sigma11 <- matrix(1:length(p1i)^2, nrow = length(p1i)) #p1ixp1i matrix
  sigma22 <- matrix(1:length(p2i)^2, nrow = length(p2i)) #p2ixp2i matrix
  sigma12 <- matrix(1:(length(p1i)*length(p2i)), nrow = length(p1i)) #p1ixp2i matrix
  
  for (i in p1i){
    for (j in p1i[column_counter:length(p1i)]){
      if (column_counter == row_counter){
        sigma11[row_counter,column_counter] <- sigma[i,j]
      }
      else{
        sigma11[row_counter,column_counter] <- sigma[i,j]
        sigma11[column_counter,row_counter] <- sigma[i,j] 
      }
      column_counter <- column_counter+1
    }
    row_counter <- row_counter + 1
    column_counter <- row_counter
  }
  
  row_counter <- 1
  column_counter <- 1
  
  for (i in p2i){
    for (j in p2i[column_counter:length(p2i)]){
      if (column_counter == row_counter){
        sigma22[row_counter,column_counter] <- sigma[i,j]
      }
      else{
        sigma22[row_counter,column_counter] <- sigma[i,j]
        sigma22[column_counter,row_counter] <- sigma[i,j] 
      }
      column_counter <- column_counter+1
    }
    row_counter <- row_counter + 1
    column_counter <- row_counter
  }
  
  row_counter <-1
  column_counter <- 1
  
  for (i in p1i){
    for (j in p2i){
      sigma12[row_counter, column_counter] <- sigma[i,j]
      column_counter <- column_counter + 1
    }
    row_counter <- row_counter + 1
    column_counter <- 1
  }
  
  sigma21 <- t(sigma12)
  
  list(mu1 = mu1, mu2 = mu2, sigma11 = sigma11, sigma22 = sigma22, sigma12 = sigma12, sigma21 = sigma21)
}

vibariate_moments <- function(x, y, joint, row_counter = 1, col_counter = 1, exy = 0){
  x_marginal <- apply(joint, 1, sum)
  y_marginal <- apply(joint, 2, sum)
  mu_x <- sum(x*x_marginal)
  mu_y <- sum(y*y_marginal)
  var_x <- sum(x_marginal*(x^2)) - mu_x^2
  var_y <- sum(y_marginal*(y^2)) - mu_y^2
  
  for (i in x){
    for (j in y){
      exy <- exy + i*j*joint[row_counter,col_counter]
      col_counter <- col_counter + 1
    }
    row_counter <- row_counter + 1
    col_counter <- 1
  }
  
  cov_ <- exy - mu_x*mu_y
  corr_ <- cov_/(sqrt(var_x*var_y))
  
  sigma_ <- matrix(c(var_x, cov_, cov_, var_y), nrow = 2)
  rho_ <- matrix(c(1, corr_, corr_, 1))
  
  list(mu_x = mu_x, mu_y = mu_y, var_x = var_x, var_y = var_y, exy = exy, cov_ = cov_, corr_ = corr_,
       sigma_ = sigma_, rho_ = rho_)
}

conditionald <- function(sigma11, sigma12, sigma22, mu1 = NULL, mu2 = NULL, x2 = NULL){
  #entrega la distribucion condicional de una particion 
  joined <- F
  variance <- sigma11 - sigma12%*% solve(sigma22) %*% t(sigma12)
  if (!is.null(mu1) && !is.null(mu2) && !is.null(x2)){
    joined <- T
    condmean <- mu1 + sigma12 %*% solve(sigma22) %*% (x2 - mu2)
    return(list(variance = variance, condmean = condmean))
  }
  mfactor <- sigma12 %*% solve(sigma22)
  list(variance = variance, m_factor = mfactor) 
}

factord <- function(mu, sigma, A){
  #entrega la distribución de una combinación lineal de va normales
  muf <- A %*% mu
  sigmaf <- A %*% sigma %*% t(A)
  list(muf = muf, sigmaf = sigmaf)
}


vectordist <- function(means, sigma, scalars){
  #distribucion de cl de vectores aleatorios
  if(!("greekLetters" %in% (.packages()))){
    library("greekLetters")
  }
  n <- dim(means[[1]])[1]
  mu <- matrix(rep(0, n), ncol = 1)
  for (i in 1:length(scalars)){
    mu <- mu + scalars[i]*means[[i]]
  }
  sigma <- sum(scalars^2)*sigma
  symbols <- paste(as.character(c(sum(scalars), sum(scalars^2))),
                   c(greeks("mu"), greeks("Sigma")), sep = "")
  list(mu = mu, sigma = sigma, symbolic = symbols, coefs = scalars)
}

#tensordist <- function(means, sigma, scalars){
  tensormean <- list()
  tensorvarcov <- list()
  for (i in 1:length(means)){
    tensormean[[i]] <- means[[i]]
  }
  names(tensormean) <- paste(rep("mu", length(means)),
                             as.character(1:length(means)), sep = "")
  counter <- 1
  for (i in 1:length(scalars)){
    for (j in 1:length(scalars)){
      sc <- sum(scalars[[i]]*scalars[[j]])
      tensorvarcov[[counter]] <- sc*sigma
      names(tensorvarcov)[counter] <- paste(i,j, sep = "")
      counter <- counter + 1
    }
  }
  list(tensormu = tensormean, tensorsigma = tensorvarcov)
#}

varitensordist <- function(vectors){
  if(!("greekLetters" %in% (.packages()))){
    library("greekLetters")
  }
  tensormean <- list()
  meanfactors <- list()
  tensorvarcov <- list()
  factors <- list()
  for (i in 1:length(vectors)){
    tensormean[[i]] <- (vectors[[i]])[[1]]
    meanfactors[[i]] <- ((vectors[[i]])[[3]])[1]
  }
  names(tensormean) <- paste(rep("mu", length(vectors)),
                             as.character(1:length(vectors)), sep = "")
  counter <- 1
  sigma <- vectors[[1]][[2]]/sum(vectors[[1]][[4]]^2)
  for (i in 1:length(vectors)){
    for (j in 1:length(vectors)){
      sc <- sum(vectors[[i]][[4]]*vectors[[j]][[4]])
      factors[[counter]] <- paste(as.character(sc), greeks("Sigma"), sep = "")
      tensorvarcov[[counter]] <- sc*sigma
      names(tensorvarcov)[counter] <- paste(i,j, sep = "")
      names(factors)[counter] <- paste(i,j, sep = "")
      counter <- counter + 1
    }
  }
  list(tensormu = tensormean, tensorsigma = tensorvarcov, 
       symbolicmu = meanfactors,symbolicsigma = factors)
}

visual <- function(tensor){
  #visualizar distribucion tensorial
  print("Vector de medias")
  for (i in 1:length(tensor[[3]])){
    print(tensor[[3]][[i]])
    print(ifelse(i == length(tensor[[3]]), "Matriz de varianzas covarianzas", "..."))
  }
  for (i in 1:length(tensor[[3]])){
    half <- tensor[[4]][[paste(i,i,sep = "")]]
    if (i > 1){
      for (k in 1:(i-1)){
        half <- paste(" ", half)
      }
    }
    if (i != length(tensor[[3]])){
      for (j in (i+1):length(tensor[[3]])){
        half <- paste(half, tensor[[4]][[paste(i,j,sep = "")]])
      }
    }
    print(half)
  }
}
