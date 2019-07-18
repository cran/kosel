cusum = function(t){
  S = cumsum(t - mean(t))
  i = which(abs(S) == max(abs(S)))
  i
}
