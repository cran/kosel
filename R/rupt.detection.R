rupt.detection = function(y,K){
  n <- length(y)
  r <- NULL
  r[K] <- n
  r[0] <- 0

  C <- matrix(nrow=n, ncol=n)
  mu <- NULL

  for(i in 0:n){
    for(j in 0:n){

      muk <- NULL
      if(i < j){
        muk <- mean(y[(i+1):j])
        C[i,j] <- sum((y[(i+1):j]- muk)^2)
      }
      else{ C[i,j] <- 999999}
    }
  }


  C2 <- matrix(nrow=(K-1), ncol=n)

  for(j in 1:(K-1)){
    for(m in 1:n){

      if(j == 1){ C2[j,m] <- C[j,m] }
      if(j != 1){
        if( m >= j){
          C2[j,m] <- 9999999

          for(h in j: m){

            if( (C2[(j-1),h] + C[h,m]) < C2[j,m]  ){
              C2[j,m] <- C2[(j-1),h] + C[h,m]
            }
          }
        }else{
          C2[j,m] <- 999999
        }
      }
    }}

  k <- K-1

  for(l in 1:k){
    TestLik <- 9999999

    for(h in (K-l): (r[K-l+1]-1)){
      if( (C2[K-l,h] + C[h,r[K-l+1]]) < TestLik ){
        TestLik <-  (C2[K-l,h] + C[h,r[K-l+1]])
        r[K-l] <- h
        mu[K+1-l] <- mean(y[(r[K-l]+1) : (r[K-l+1]) ] )
      }
    }
  }

  mu[1] <- mean(y[1:r[1]])
  r <- r[-K]


  return(r)
}
