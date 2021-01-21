# A set of functions to create constraint matrices.

# Create a finite difference matrix, any order d including fractional and negative
FiniteDiffM = function(N,d){
  # Find finite difference coefficient
  u=1:(N-1)
  wt = rev(c(1,cumprod ((u-d-1)/u))) # generalized binomial series
  Dmat = Matrix::Matrix(0,N,N,sparse=TRUE)
  # shift through, create matrix
  for(rw in 1:N){
    Dmat[rw,1:rw] = wt[(N-rw+1):N]
  }

  # clear the first rows
  if(d>0){
    D = ceiling(d);
    Dmat[1:D,1:D] = 0
  }
  return(Dmat)
}

# Create an Exponential Smoothing Constraint
ExpSmMat = function(N,alpha,double=FALSE){
  alpha = min(alpha,0.99999)
  alpha = max(alpha,0)
  k=N:1
  wt=(1-alpha)^k
  mat = matrix(0,N,N)
  mat[1,1]=1
  for(rw in 2:N){
    mat[rw,1:(rw-1)] = wt[(N-rw+2):N]
  }
  mat = mat/rowSums(mat)

  if(double){
    mat = 2*mat-mat%*%mat # double (Browns) exponential smoothing
    mat[2,1:2] = c(0,1) # don't constrain initial slope
  }

  # put in difference form
  ESM = diag(N)-mat
  return(ESM)
}


# AR matrix
ARmat = function(N,ar){

  # initial stuff
  mat = diag(N) # 1 - lagged...
  lngAR = length(ar)
  indset = (1:lngAR)-lngAR-1

  # Not possible in this case
  if(lngAR>=N){
    return(0*mat)
  }

  # cascade down
  for(rw in (lngAR+1):N){
    mat[rw,indset+rw] = -ar
  }
  return(mat)
}



# Create seasonal difference matrix
Seasonal_DM= function(N,lag=12,sumFirst=TRUE){
  Dm = diag(N)
  if( N <= lag){
    return(0*Dm)
  }
  Dm[row(Dm)==(col(Dm)+lag)]=-1
  Dm[1:lag,1:lag]=0
  if(sumFirst){
    Dm[1,1:lag] =1
  }
  return(Dm)
}


