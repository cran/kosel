#' Variable selection with the knockoffs procedure.
#'
#' @description Performs variable selection from an object (vector of statistics W) returned by \code{\link{ko.glm}} or \code{\link{ko.ordinal}}.
#'
#' @param W A vector of length nvars corresponding to the statistics W. Object returned by the functions \code{ko.glm} or \code{ko.ordinal}.
#' @param print Logical. If \code{TRUE}, positive statistics W are displayed in increasing order. If \code{FALSE}, nothing is displayed. If \code{method = 'manual'}, \code{print} is automatically \code{TRUE}.
#' @param method Can be \code{'stats'}, \code{'gaps'} or \code{'manual'}. If \code{'stats'}, the threshold used is the W-threshold. If \code{'gaps'}, the threshold used is the gaps-threshold. If \code{'manual'}, the user can choose its own threshold using the graph of the positive statistics W sorted in increasing order.
#'
#' @return A list containing two elements:
#' \itemize{
#' \item \code{threshold}  A positive real value corresponding to the threshold used.
#' \item \code{estimation} A binary vector of length nvars corresponding to the variable selection: 1*(W >= threshold). 1 indicates that the associated covariate belongs to the estimated model.
#' }
#'
#' @seealso \code{\link{ko.glm}}, \code{\link{ko.ordinal}}
#'
#' @export ko.sel
#' @importFrom graphics abline plot text
#'
#' @references Gegout-Petit Anne, Gueudin Aurelie, Karmann Clemence (2019). \href{https://arxiv.org/pdf/1907.03153.pdf}{\emph{The revisited knockoffs method for variable selection in L1-penalised regressions}, arXiv:1907.03153.}
#'
#'
#' @examples
#'
#' library(graphics)
#'
#' # linear Gaussian regression
#' n = 100
#' p = 20
#' set.seed(11)
#' x = matrix(rnorm(n*p),nrow = n,ncol = p)
#' beta = c(rep(1,5),rep(0,15))
#' y = x%*%beta + rnorm(n)
#' W = ko.glm(x,y)
#' ko.sel(W, print = TRUE)
#'
#'
#' # logistic regression
#' n = 100
#' p = 20
#' set.seed(11)
#' x = matrix(runif(n*p, -1,1),nrow = n,ncol = p)
#' u = runif(n)
#' beta = c(c(3:1),rep(0,17))
#' y = rep(0, n)
#' a = 1/(1+exp(0.1-x%*%beta))
#' y = 1*(u>a)
#' W = ko.glm(x,y, family = 'binomial', nVal = 50)
#' ko.sel(W, print = TRUE)
#'
#'
#' # cumulative logit regression
#' n = 100
#' p = 10
#' set.seed(11)
#' x = matrix(runif(n*p),nrow = n,ncol = p)
#' u = runif(n)
#' beta = c(3,rep(0,9))
#' y = rep(0, n)
#' a = 1/(1+exp(0.8-x%*%beta))
#' b = 1/(1+exp(-0.6-x%*%beta))
#' y = 1*(u<a) + 2*((u>=a) & (u<b)) + 3*(u>=b)
#' W = ko.ordinal(x,as.factor(y), nVal = 20)
#' ko.sel(W, print = TRUE)
#'
#'
#' # adjacent logit regression
#' n = 100
#' p = 10
#' set.seed(11)
#' x = matrix(rnorm(n*p),nrow = n,ncol = p)
#' U = runif(n)
#' beta = c(5,rep(0,9))
#' alpha = c(-2,1.5)
#' M = 2
#' y = rep(0, n)
#' for(i in 1:n){
#'   eta = alpha + sum(beta*x[i,])
#'   u = U[i]
#'   Prob = rep(1,M+1)
#'   for(j in 1:M){
#'    Prob[j] = exp(sum(eta[j:M]))
#'   }
#'   Prob = Prob/sum(Prob)
#'   C = cumsum(Prob)
#'   C = c(0,C)
#'   j = 1
#'   while((C[j]> u) || (u >= C[j+1])){j = j+1}
#'   y[i] = j
#' }
#' W = ko.ordinal(x,as.factor(y), family = 'acat', nVal = 10)
#' ko.sel(W, method = 'manual')
#' 0.4
#'
#'
#' # How to use randomness?
#' n = 100
#' p = 20
#' set.seed(11)
#' x = matrix(rnorm(n*p),nrow = n,ncol = p)
#' beta = c(5:1,rep(0,15))
#' y = x%*%beta + rnorm(n)
#' Esti = 0
#' for(i in 1:100){
#'   W = ko.glm(x,y, random = TRUE)
#'   Esti = Esti + ko.sel(W, method = 'gaps')$estimation
#' }
#' Esti
#'
ko.sel = function(W, print = FALSE, method = 'stats'){
  Test = is.matrix(W)*(nrow(W) == 1)
  if(!is.vector(W) && !Test){stop("W should be a vector of length nvars.")}
  if(method != 'stats' && method != 'gaps' && method != 'manual'){stop("method should be 'stats', 'gaps' or 'manual.")}
  if(print != FALSE && print != TRUE){stop("print should be TRUE or FALSE.")}
  #########

  p = length(W)
  WW = sort(W[W > 0])
  if(method == 'stats'){
    tt = min(rupt.detection(WW,2),cusum(WW))
    if(print == TRUE){
      cc = NULL
      U = unique(WW)
      for(i in 1:length(U)){
        cc = c(cc, which(W == U[i]))
      }
      plot(WW, type = 'n', xlab = '', ylab = 'Positive statistics W')
      text(1:length(WW),WW,cc, cex = 0.7)
      abline(h = WW[tt+1], col = 'red')
    }
    t = WW[tt+1]
    Esti= 1*(W > WW[tt])
  }
  if(method == 'gaps'){
    l = length(WW)
    if(l == 0){Esti = rep(0,p)}
    if(l == 1){Esti = 1*(W > 0)}
    if(l > 1){
      e = rep(0, l - 1)
      for (i in 1:length(e)){
        e[i] = WW[i+1] - WW[i]
      }
      tt = min(rupt.detection(e,2), cusum(e), cusum(-e), rupt.detection(-e,2))
      if(print == TRUE){
        cc = NULL
        U = unique(WW)
        for(i in 1:length(U)){
          cc = c(cc, which(W == U[i]))
        }
        plot(WW, type = 'n', xlab = '', ylab = 'Positive statistics W')
        text(1:length(WW),WW,cc, cex = 0.7)
        abline(h = WW[tt+2], col = 'red')
      }
      t = WW[tt+2]
      Esti= 1*(W > WW[tt+1])
    }
  }
  if(method == 'manual'){
    cc = NULL
    U = unique(WW)
    for(i in 1:length(U)){
      cc = c(cc, which(W == U[i]))
    }
    plot(WW, type = 'n', xlab = '', ylab = 'Positive statistics W')
    text(1:length(WW),WW,cc, cex = 0.7)
    print('Enter a positive real number as a threshold:')
    t = scan(what = numeric(), nmax = 1,nlines = 1)
    abline(h = t, col = 'red')
    Esti = 1*(W >= t)
  }
  list(estimation = Esti, threshold = t)
}
