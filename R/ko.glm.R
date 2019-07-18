#' Statistics of the knockoffs procedure for glmnet regression models.
#'
#' @description Returns the vector of statistics W of the revisited knockoffs procedure for regressions available in the R package \code{glmnet}. Most of the parameters come from \code{glmnet()}. See \href{https://CRAN.R-project.org/package=glmnet}{\code{glmnet} documentation} for more details.
#'
#' @param x Input matrix, of dimension nobs x nvars; each row is an observation vector. Can be in sparse matrix format (inherit from class "\code{sparseMatrix}" as in package \code{Matrix}; not yet available for \code{family="cox"})
#' @param y Response variable. Quantitative for \code{family="gaussian"}, or \code{family="poisson"} (non-negative counts). For family="binomial" should be either a factor with two levels, or a two-column matrix of counts or proportions (the second column is treated as the target class; for a factor, the last level in alphabetical order is the target class). For \code{family="multinomial"}, can be a \code{nc>=2} level factor, or a matrix with \code{nc} columns of counts or proportions. For either \code{"binomial"} or \code{"multinomial"}, if \code{y} is presented as a vector, it will be coerced into a factor. For \code{family="cox"}, \code{y} should be a two-column matrix with columns named 'time' and 'status'. The latter is a binary variable, with '1' indicating death, and '0' indicating right censored. The function \code{Surv()} in package survival produces such a matrix.
#' @param family Response type: "gaussian","binomial","poisson","multinomial","cox". Not available for "mgaussian".
#' @param alpha The elasticnet mixing parameter, with 0 <= \code{alpha} <= 1. \code{alpha=1} is the lasso penalty, and \code{alpha=0} the ridge penalty. The default is 1.
#' @param type.gaussian See \code{glmnet} documentation.
#' @param type.logistic See \code{glmnet} documentation.
#' @param type.multinomial See \code{glmnet} documentation.
#' @param nVal Length of lambda sequence - default is 50.
#' @param random If \code{TRUE}, the matrix of knockoffs is different for every run. If \code{FALSE}, a seed is used so that the knockoffs are the same. The default is \code{FALSE}.
#'
#' @return A vector of dimension nvars corresponding to the statistics W.
#'
#' @seealso \code{\link{ko.sel}}
#' @export ko.glm
#' @importFrom glmnet glmnet
#'
#' @examples
#' # see ko.sel
#'
ko.glm = function(x,y,family = "gaussian", alpha = 1, type.gaussian = ifelse(nvars<500,"covariance","naive"), type.logistic = "Newton", type.multinomial = "ungrouped", nVal = 50, random = FALSE)
{
  # Bonne utilisation :
  if (nVal < 2){stop("nVal should be an integer >= 2.")}
  if(family == 'mgaussian'){stop("Variable selection not available for family 'mgaussian'.")}
  if(random != FALSE && random != TRUE){stop("random should be TRUE or FALSE.")}
  nVal = ceiling(nVal)
  ####################

  p = ncol(x)
  n = nrow(x)
  nvars = n
  K = length(unique(y))
  if(random == FALSE){set.seed(1)}
  KO = x[sample(1:n, n),]
  L = glmnet(cbind(x,KO),y, family = family, alpha = alpha, type.gaussian = type.gaussian, type.logistic = type.logistic, type.multinomial = type.multinomial)$lambda
  m = min(L)
  M = max(L)
  Lambda = seq(M,m, length = nVal)
  G = glmnet(cbind(x,KO),y,family = family, alpha = alpha, lambda = Lambda, type.gaussian = type.gaussian, type.logistic = type.logistic, type.multinomial = type.multinomial)$beta
  if(family == 'multinomial'){
    K = length(unique(y))
    Z = 0
    G = lapply(G, function(x){t(1*(as.matrix(x) != 0))})
    for(i in 1:K){
      Z = Z + G[[i]]
    }
  }
  else{
    Z = as.matrix(G)
    Z = t(Z)
    Z = 1*(Z != 0)
  }
  Z = apply(Z,2,function(vec){vec*Lambda})
  Z = apply(Z,2,max)
  Z = t(Z)

  colnames(Z) = 1:(2*p)
  W = rep(0,p)
  for (i in 1:p){
    if (Z[i] > Z[i+p]){
      W[i] = Z[i]
    }
    if (Z[i] <= Z[i+p])
    {
      W[i] = -Z[i+p]
    }
  }
  W = t(W)
  colnames(W) = 1:p

  W
}
