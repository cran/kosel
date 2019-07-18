#' Statistics of the knockoffs procedure for ordinalNet regression models.
#'
#' @description Returns the vector of statistics W of the revisited knockoffs procedure for regressions available in the R package \code{ordinalNet}. Most of the parameters come from \code{ordinalNet()}. See \href{https://CRAN.R-project.org/package=ordinalNet}{\code{ordinalNet} documentation} for more details.
#'
#' @param x Covariate matrix, of dimension nobs x nvars; each row is an observation vector. It is recommended that categorical covariates are converted to a set of indicator variables with a variable for each category (i.e. no baseline category); otherwise the choice of baseline category will affect the model fit.
#' @param y Response variable. Can be a factor, ordered factor, or a matrix where each row is a multinomial vector of counts. A weighted fit can be obtained using the matrix option, since the row sums are essentially observation weights. Non-integer matrix entries are allowed.
#' @param family Specifies the type of model family. Options are "cumulative" for cumulative probability, "sratio" for stopping ratio, "cratio" for continuation ratio, and "acat" for adjacent category.
#' @param reverse Logical. If TRUE, then the "backward" form of the model is fit, i.e. the model is defined with response categories in reverse order. For example, the reverse cumulative model with K+1 response categories applies the link function to the cumulative probabilities P(Y >= 2), …, P(Y >= K+1), rather then P(Y <= 1), …, P(Y <= K).
#' @param link Specifies the link function. The options supported are logit, probit, complementary log-log, and cauchit.
#' @param alpha The elastic net mixing parameter, with \code{0 <= alpha <= 1}. \code{alpha=1} corresponds to the lasso penalty, and \code{alpha=0} corresponds to the ridge penalty.
#' @param parallelTerms Logical. If \code{TRUE}, then parallel coefficient terms will be included in the model. \code{parallelTerms} and \code{nonparallelTerms} cannot both be \code{FALSE}.
#' @param nonparallelTerms Logical. if \code{TRUE}, then nonparallel coefficient terms will be included in the model. \code{parallelTerms} and \code{nonparallelTerms} cannot both be \code{FALSE}. Default is \code{FALSE}. \code{nonparallelTerms = TRUE} is highly discouraged.
#' @param nVal Length of lambda sequence - default is 100.
#' @param warn Logical. If \code{TRUE}, the following warning message is displayed when fitting a cumulative probability model with \code{nonparallelTerms=TRUE} (i.e. nonparallel or semi-parallel model). "Warning message: For out-of-sample data, the cumulative probability model with \code{nonparallelTerms=TRUE} may predict cumulative probabilities that are not monotone increasing." The warning is displayed by default, but the user may wish to disable it.
#' @param random If \code{TRUE}, the matrix of knockoffs is different for every run. If \code{FALSE}, a seed is used so that the knockoffs are the same. The default is \code{FALSE}.
#'
#' @note \code{nonparallelTerms = TRUE} is highly discouraged because the knockoffs procedure does not suit well to this setting.
#'
#' @return A vector of dimension nvars corresponding to the statistics W.
#' @seealso \code{\link{ko.sel}}
#' @export ko.ordinal
#' @importFrom ordinalNet ordinalNet
#'
#' @examples
#' # see ko.sel
#'
#'
ko.ordinal = function(x,y,family = 'cumulative', reverse = FALSE, link = 'logit', alpha = 1, parallelTerms = TRUE, nonparallelTerms = FALSE, nVal = 100, warn = FALSE, random = FALSE)
{
  # Bonne utilisation :
  if (!is.matrix(x)){stop("x should be a matrix.")}
  if (nVal < 2){stop("nVal should be an integer >= 2.")}
  if(random != FALSE && random != TRUE){stop("random should be TRUE or FALSE.")}
  nVal = ceiling(nVal)
  ####################

  p = ncol(x)
  n = nrow(x)
  K = length(unique(y))
  if(random == FALSE){set.seed(1)}
  KO = x[sample(1:n, n),]
  L = ordinalNet(cbind(x,KO),y, reverse = reverse, family = family, link = link, alpha = alpha, parallelTerms = parallelTerms, nonparallelTerms = nonparallelTerms, warn = warn)$lambdaVals
  m = min(L)
  M = max(L)
  Lambda = seq(M,m, length = nVal)
  coef = ordinalNet(cbind(x,KO),y,family = family, reverse = reverse, link = link, alpha = alpha, parallelTerms = parallelTerms, nonparallelTerms = nonparallelTerms,lambdaVals = Lambda, warn = warn)$coefs
  l = ncol(coef)
  Z = coef[,(-(1:(K-1)))]

  # Calcul des statistiques
  Z = 1*(Z != 0)
  if(parallelTerms == TRUE && nonparallelTerms == FALSE){
    Z = apply(Z,2,function(vec){vec*Lambda})
    Z = apply(Z,2,max)
    Z = t(Z)
  }
  if(parallelTerms == FALSE && nonparallelTerms == TRUE){
    Z = apply(Z, 1, function(x){f = 0
    for(i in 1:(K-1)){
      f = f + x[((i-1)*2*p+1):(i*2*p)]
    }
    f})
    Z = t(Z)
    Z = 1*(Z != 0)
    Z = apply(Z,2,function(vec){vec*Lambda})
    Z = apply(Z,2,max)
    Z = t(Z)
  }
  if(parallelTerms == TRUE && nonparallelTerms == TRUE){
    Z = apply(Z, 1, function(x){f = 0
    for(i in 1:K){
      f = f + x[((i-1)*2*p+1):(i*2*p)]
    }
    f})
    Z = t(Z)
    Z = 1*(Z != 0)
    Z = apply(Z,2,function(vec){vec*Lambda})
    Z = apply(Z,2,max)
    Z = t(Z)
  }
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
