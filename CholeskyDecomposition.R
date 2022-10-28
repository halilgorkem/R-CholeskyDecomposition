A <- matrix(c(4,-1,1,
              -1,4.25,2.75,
              1,2.75,3.50), ncol = 3, nrow = 3, byrow = T)
cholesky <- function(A)
{
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
  list(L=L)
}
cholesky(A)
