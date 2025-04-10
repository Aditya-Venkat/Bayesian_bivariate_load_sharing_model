#Simulate Data
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


#Function to find MLEs of lambda given alpha
find_MLE_lambdas <- function(dat , alpha){
  
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
  
  
  
  lambda1 <- nrow(set1)/(alpha*(sum(set3[,1]) + sum(set1[,1]) + sum(set2[,2])))
  lambda2 <- nrow(set2)/(alpha*(sum(set3[,1]) + sum(set1[,1]) + sum(set2[,2])))
  lambda3 <- nrow(set3)/(alpha*(sum(set3[,1]) + sum(set1[,1]) + sum(set2[,2])))
  lambda4 <- nrow(set2)/(alpha*(sum(set2[,1]) - sum(set2[,2])))
  lambda5 <- nrow(set1)/(alpha*(sum(set1[,2]) - sum(set1[,1])))
  
  return(c(lambda1, lambda2, lambda3, lambda4, lambda5))
}

#1D optimization to find alpha, fails since the function is strictly decreasing
find_MLE_alpha <- function(dat){
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
  
  l2 <- function(alpha){
    return((n3 + 2*n1 + 2*n2)*log(alpha))
  }
  
  h_alpha <- function(alpha){
    lambdas <- find_MLE_lambdas(dat, alpha)
    h <- n3*log(lambdas[3]) + n1*log(lambdas[1]) + n1*log(lambdas[5]) + n2*log(lambdas[2]) + n2*log(lambdas[4]) -
      (lambdas[1] + lambdas[2] + lambdas[3])*alpha*sum(set3[,1]) - 
      (lambdas[1] + lambdas[2] + lambdas[3])*alpha*sum(set1[,1]) +lambdas[5]*(-alpha*sum(set1[,2]) + alpha*sum(set1[,1])) -
      (lambdas[1] + lambdas[2] + lambdas[3])*alpha*sum(set2[,2]) + lambdas[4]*(-alpha*sum(set2[,1]) + alpha*sum(set2[,2]))
      + l2(alpha)
    return(h)
  }
  
  result <- optim(1 , h_alpha , method = "Brent" ,lower = 0, upper= 10)
  
 
  return(result$par) 
}


#Computation of MLEs of lambdas as errors
final <- matrix(nrow = 0 , ncol = 5)

for(n in c(50, 1e2,1e3)){
  
  estims <- matrix(nrow = 0, ncol = 5)
  errors <- matrix(nrow = 0 , ncol = 5)
  for(i in 1:1000){
    print(i)
    dat <- draw_biv_data(n, c(5,6,7,6,5) , 2)
    lambs <- find_MLE_lambdas(dat,2)
    estims <- rbind(estims, lambs)
    errors <- rbind(errors, (c(5,6,7,6,5) - lambs)^2)
  }
  estims <- colMeans(estims)
  errors <- colMeans(errors)
  final <- rbind(final , estims)
  final <- rbind(final , errors)
  
}


