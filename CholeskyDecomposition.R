A <- matrix(c(4,-1,1,
              -1,4.25,2.75,
              1,2.75,3.50), ncol = 3, nrow = 3, byrow = T)
b <- c(1,1,1)
cholesky_decomposition <- function(A,b)
{
  X <- rep(0,n)
  Y <- rep(0,n)
  n <- nrow(A)
  L <- matrix(c(0), n,n)
  
  #step1
  L[1,1] <- sqrt(A[1,1])
  L
  #step2
  for(j in 2:n)
  {
    L[j,1] <- A[j,1]/L[1,1]
  }
  L
  #step3
  for(i in 2:(n-1))
  {
    #step4
    L[i,i] <- sqrt(A[i,i] - sum((L[i, 1:(i-1)])^2))
    L
    #step5
    for(j in (i+1):n)
    {
      L[j,i] = (A[j,i] - sum(L[j,1:(i-1)]*L[i,1:(i-1)]))/L[i,i]
    }
  }
  L
  #step6
  L[n,n] <- sqrt(A[n,n] - sum((L[n, 1:(n-1)])^2))
  
  #step7
  #list(L=L)
  
  #step8
  Y <- b[1]/L[1,1]
  
  #step9
  for(i in 2:n)
  {
    Y[i] <- (b[i] - sum(L[i,1:(i-1)]*Y[1:(i-1)]))/L[i,i]
  }
  #step10
  X[n] <- Y[n]/L[n,n]
  
  #step11
  for(i in (n-1):1)
  {
    X[i] <- (Y[i] - sum(L[(i+1):n,i]*X[(i+1):n]))/L[i,i]
  }
  #step12
  list(X=X,L=L)
}
cholesky_decomposition(A,b)
