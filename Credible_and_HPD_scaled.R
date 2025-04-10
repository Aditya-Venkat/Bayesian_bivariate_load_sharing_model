library(matrixStats)

#simulate data
draw_biv_data <- function(n,lambdas, alpha ){
  
  out <- matrix(nrow = n , ncol = 2)
  u1 <- runif(n , 0 , 1)
  u2 <- runif(n, 0, 1)
  u3 <- runif(n, 0, 1)
  
  v1 <- -log(u1^(1/lambdas[1]))/alpha
  v2 <- -log(u2^(1/lambdas[2]))/alpha
  v3 <- -log(u3^(1/lambdas[3]))/alpha
  
  x <- sapply(1:n, function(i){min(v1[i],v3[i])})
  y <- sapply(1:n, function(i){min(v2[i],v3[i])})
  
  final_draw <- function(v3,x,y,lambdas){
    if(x  == y){
      return(c(v3,v3))
    }
    else if(x > y){
      y2 <- y
      u <- runif(1,0,1)
      y1 <- -log( (u^(1/lambdas[4]))*(exp(-alpha*y2)) )/alpha
      return(c(y1,y2))
    }
    else{
      y1 <- x
      u <- runif(1,0,1)
      y2 <- -log( (u^(1/lambdas[5]))*(exp(-alpha*y1)) )/alpha
      return(c(y1,y2))
      
    }
  }
  
  out <- sapply(1:n ,function(i){final_draw(v3[i],x[i],y[i],lambdas)})
  return(t(out))
}

## assuming priors of alpha as gamma(a,b) and of lambdas as gamma(ci,di)
#Importance function (scaled)
post_lambda_samp <- function(n,dat, a, b, cs, ds){
  
  for(i in 1:5) {
    assign(paste0("c", i), cs[i], envir = environment())
    assign(paste0("d", i), ds[i], envir = environment())
  }
  
  set1 <- matrix(nrow = 0, ncol = 2)
  set2 <- matrix(nrow = 0, ncol = 2)
  set3 <- matrix(nrow = 0, ncol = 2)
  
  
  split_data <- apply(dat, 1, function(row) {
    if (row[1] < row[2]) {
      set1 <<- rbind(set1, (row))
    } else if (row[1] > row[2]) {
      set2 <<- rbind(set2, (row))
    } else {
      set3 <<- rbind(set3, (row))
    }
  })
  n1 <- nrow(set1)
  n2 <- nrow(set2)
  n3 <- nrow(set3)
  
  thetas <- c((n2+c1)/d1 , (n3 + c2)/d2 , (n1 + c3)/d3 , (n3 + c4)/d4 , (n2 + c5)/d5)
  
  lambda1 <- rgamma(n, shape = n2 + c1 , rate = d1)
  lambda2 <- rgamma(n , shape = n3 + c2 , rate = d2)
  lambda3 <- rgamma(n , shape = n1 + c3 , rate = d3)
  lambda4 <- rgamma(n , shape = n3 + c4 , rate = d4)
  lambda5 <- rgamma(n , shape = n2 + c5 , rate = d5)
  
  lambdas <- cbind(lambda1/thetas[1], lambda2/thetas[2], lambda3/thetas[3], lambda4/thetas[4], lambda5/thetas[5])
  
  return(lambdas)
  
  
}

#Code to compute credible intervals of lambdai
post_credible <- function(i,alpha,dat, a, b, cs, ds, n = 1e3){
  
  set1 <- matrix(nrow = 0, ncol = 2)
  set2 <- matrix(nrow = 0, ncol = 2)
  set3 <- matrix(nrow = 0, ncol = 2)
  
  
  split_data <- apply(dat, 1, function(row) {
    if (row[1] < row[2]) {
      set1 <<- rbind(set1, (row))
    } else if (row[1] > row[2]) {
      set2 <<- rbind(set2, (row))
    } else {
      set3 <<- rbind(set3, (row))
    }
  })
  n1 <- nrow(set1)
  n2 <- nrow(set2)
  n3 <- nrow(set3)
  
  # if(is.null(ds) || is.null(cs)){
  #   print("Setting values of cs and ds")
  #   ds <- c(n2,n3,n1,n3,n2)
  #   cs <- ds*sqrt(ds)
  #   
  # }
  
  for(j in 1:5){
    assign(paste0("c", j), cs[j], envir = environment())
    assign(paste0("d", j), ds[j], envir = environment())
  }
  
  thetas <- c((n2+c1)/d1 , (n3 + c2)/d2 , (n1 + c3)/d3 , (n3 + c4)/d4 , (n2 + c5)/d5)
  
  lambdas <- post_lambda_samp(n , dat, a, b , cs, ds)
  
  lambdas <- lambdas[order(lambdas[,i]),]
  
  
  S <- (lambdas[,1] + lambdas[,2] + lambdas[,3])*sum(set1[,1]) + 
    (lambdas[,1] + lambdas[,2] + lambdas[,3] - lambdas[,5])*sum(set2[,1]) + lambdas[,5]*sum(set2[,2]) +
    (lambdas[,1] + lambdas[,2] + lambdas[,3] - lambdas[,4])*sum(set3[,2]) + lambdas[,4]*sum(set3[,1])
  
  log_w <- lgamma(n2+c1) + lgamma(n3+c2) + lgamma(n1+c3) + lgamma(n3+c4) + lgamma(n2+c5)-
    ((n2+c1)*log(d1) + (n3+c2)*log(d2) + (n1+c3)*log(d3) + (n3+c4)*log(d4) + (n2+c5)*log(d5) +
       (a + n1 + 2*n2 + 2*n3)*log(b+S) + (n2 + c1 )*log(thetas[1]) + (n3 + c2 )*log(thetas[2])
     + (n1 + c3 )*log(thetas[3]) + (n3 + c4 )*log(thetas[4]) + (n2 + c5 )*log(thetas[5])) - ((1-thetas[1])*d1*lambdas[1] + (1 - thetas[2])*d2*lambdas[2] + (1-thetas[3])*d3*lambdas[3] + 
                                                                                                       (1-thetas[4])*d4*lambdas[4] + (1-thetas[5])*d5*lambdas[5])
  logsumw <- logSumExp(log_w)
  
  
  # w <- (gamma(n2+c1) * gamma(n3+c2) * gamma(n1+c3) * gamma(n3 + c4) * gamma(n2+c5))/
  #   (d1^(n2+c1) * d2^(n3+c2) * d3^(n1 + c3) * d4^(n3+c4) * d5^(n2+c5) *(b + S)^(a + n1 + 2*n2 + 2*n3) )
  w <- exp(log_w - logsumw)
  
  find_theta_alpha <- function(alpha){
    cum_w <- cumsum(w)
    if(alpha == 0){
      return(lambdas[1,i])
    }
    else{
      idx <- which(alpha <= cum_w)
      ind <- idx[1]
      return(lambdas[ind,i])
    }
  }
  
  
  return(c(find_theta_alpha(alpha/2) , find_theta_alpha(1-alpha/2)))
  
}

#Code to compute HPD intervals of lambdai
post_HPD <- function(i,alpha,dat, a, b, cs = NULL, ds = NULL , n = 50){
  
  ## CODE A CHECK FOR UNIMODALITY CHECKING
  set1 <- matrix(nrow = 0, ncol = 2)
  set2 <- matrix(nrow = 0, ncol = 2)
  set3 <- matrix(nrow = 0, ncol = 2)
  
  
  split_data <- apply(dat, 1, function(row) {
    if (row[1] < row[2]) {
      set1 <<- rbind(set1, (row))
    } else if (row[1] > row[2]) {
      set2 <<- rbind(set2, (row))
    } else {
      set3 <<- rbind(set3, (row))
    }
  })
  n1 <- nrow(set1)
  n2 <- nrow(set2)
  n3 <- nrow(set3)
  
  if(is.null(ds) || is.null(cs)){
    print("Setting values of cs and ds")
    ds <- c(n2,n3,n1,n3,n2)
    cs <- ds*sqrt(ds)
    
  }
  
  
  for(j in 1:5) {
    assign(paste0("c", j), cs[j], envir = environment())
    assign(paste0("d", j), ds[j], envir = environment())
  }
  
  
  n <- dim(dat)[1]
  CIS <- matrix(nrow = 0 , ncol = 2)
  
  find_theta_alpha <- function(alpha){
    cum_w <- cumsum(w)
    if(alpha == 0){
      return(lambdas[1,i])
    }
    else{
      idx <- which(alpha <= cum_w)
      ind <- idx[1]
      return(lambdas[ind,i])
    }
  }
  for(j in 1:(n - floor((1-alpha)*n))){
    
    CI <- c(find_theta_alpha(j/n) , find_theta_alpha((j + floor((1-alpha)*n))/n))
    CIS <- rbind(CIS,CI)
    
  }
  
  widths <- CIS[,2] - CIS[,1]
  ind <- which.min(widths)
  
  HBD <- CIS[ind,]
  
  return(HBD)
  
}


dat <- draw_biv_data(50, c(5,6,7,6,5) , 7)
post_credible(1, 0.05, dat, 2,0.3, c(2,3,4,3,2) , c(0.5,0.5,0.5,0.5,0.5))
