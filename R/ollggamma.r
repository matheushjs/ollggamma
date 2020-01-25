
#' Generalized Gamma Probability Distribution
#'
#' Fast implementation of density, distribution function, quantile function
#' and random generation for the Generalized Gamma probability distribution.
#'
#' @param x,q	          vector of quantiles.
#' @param p	              vector of probabilities.
#' @param n	              number of observations. If \code{length(n) > 1},
#'                        the length is taken to be the number required.
#' @param a,b,k	          Parameters of the distribution, all of which must be positive.
#' @param log,log.p	      logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail	  logical; if TRUE (default), probabilities are \eqn{P[X \le x]}
#'                        otherwise, \eqn{P[X > x]}.
#'
#' @details
#'
#' The generalized gamma distribution proposed by Stacy (1962) has parameters
#' \eqn{a, d, p}, but here we adopt the reparametrization
#' \deqn{
#'   a = a
#' }
#' \deqn{
#'   b = p
#' }
#' \deqn{
#'   k = \frac{d}{p}
#' }
#' as is used by the R package *flexsurv*.
#'
#' Probability density function
#' \deqn{
#'    f(x) = \frac{b x^{bk-1} \exp[-(x/a)^b]}{a^{bk} \Gamma(k)}
#' }
#'
#' Cumulative density function
#' \deqn{
#'    F(x) = \frac{\gamma(k, (x/a)^b)}{\Gamma(k)}
#' }
#'
#' The above function can be written in terms of a \eqn{Gamma(\alpha, \beta)}.
#' Let \eqn{T \sim Gamma(k, 1)} and its cumulative distribution be denoted as \eqn{F_T(t)},
#' then the cumulative density function of the generalized gamma distribution can be
#' written as
#' \deqn{
#'    F(x) = F_T( (x/a)^b )
#' }
#' which allows us to write the quantile function of the generalized gamma in terms of
#' the gamma one (\eqn{Q_T(u)} is the quantile function of \eqn{T})
#' \deqn{
#'    Q(u) = (Q_T(u) \cdot a)^{1/b}
#' }
#' from which random numbers can be drawn.
#'
#' @references
#' Stacy, E. W. (1962). A generalization of the gamma distribution.
#' The Annals of mathematical statistics, 33(3), 1187-1192.
#'
#' @name G.Gamma
#' @aliases Generalized-Gamma
#' @aliases GGamma
#'
#' @keywords distribution
#' @keywords univar
#' @keywords models
#' @keywords survival
#' @concept Univariate
#' @concept Continuous
#' @concept Lifetime
#'
#' @examples
#'
#' x = seq(0.001, 5, length=1000);
#' plot(x, dggamma(x, 3, 1.8, 0.5), col=2, type="l", lwd=4, ylim=c(0, 1));
#' lines(x, pggamma(x, 3, 1.8, 0.5), col=4, type="l", lwd=4, ylim=c(0, 1));
#' legend("right", c("PDF", "CDF"), col=c(2, 4), lwd=4);
#'
#' r = rgamma(n = 100, 2, 2);
#' lik = function(params) -sum(dggamma(r, params[1], params[2], params[3], log=TRUE));
#' optPar = optim(lik, par=c(1, 1, 1), method="L-BFGS", lower=0.00001, upper=Inf)$par;
#' x = seq(0.001, 5, length=1000);
#' plot(x, dgamma(x, 2, 2), type="l", col=2, lwd=4, ylim=c(0, 1));
#' lines(x, dggamma(x, optPar[1], optPar[2], optPar[3]), col=4, lwd=4);
#' legend("topright", c("Gamma(shape=2, rate=2)", "MLE Gen. Gamma"), col=c(2, 4), lwd=4);
#'
#' @export

dollggamma = function(x, a, b, k, lambda, log=F){
	cdf  = pgamma((x / a)**b, shape=k, rate=1);                             # This is pggamma
	lpdf = log(b) - lgamma(k) + (b*k - 1)*log(x) - (b*k)*log(a) - (x/a)**b; # This is log(dggamma)
	result = log(lambda) + lpdf + (lambda - 1)*(log(cdf) + log(1 - cdf)) - 2*log(cdf**lambda + (1 - cdf)**lambda);
	result[cdf == 1] = log(0); # We analyzed the formulas and concluded this is right.
	if(!log) result = exp(result);
	return(result);
}

#' @rdname G.Gamma
#' @export

pollggamma = function(q, a, b, k, lambda, lower.tail = TRUE, log.p = FALSE){
	cdf = pgamma( (q / a)**b, shape=k, rate=1); # This is pggamma
	cdf = 1 / (1 + (1/cdf - 1)**lambda)     # now apply the odd-logistic
	if(!lower.tail) cdf = 1 - cdf;
	if(log.p) cdf = log(cdf);
	return(cdf);
}

#' @rdname G.Gamma
#' @export

qollggamma = function(p, a, b, k, lambda, lower.tail = TRUE, log.p = FALSE){
	if(log.p) p = exp(p);
	if(!lower.tail) p = 1 - p;
	quantile = 1 / ((1/p - 1)**(1/lambda) + 1)
	quantile = qgamma(p, shape=k, rate=1); # This is qggamma
	quantile = a * quantile**(1/b);        # This is also qggamma
	return(quantile);
}


#' @rdname G.Gamma
#' @export

rollggamma = function(n, a, b, k, lambda){
	return( qollggamma(runif(0, 1, n=n), a, b, k, lambda) );
}
