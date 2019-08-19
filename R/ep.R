
# functions -------------------------------------------------------------------
tRoot <- function(delta, test_stat, df, conf_level) {
  #' This function is used to find the root for a t-distribution pivotal quantity 
  #'
  #' This function returns the difference between the lower tail
  #' probability of a non-central t-distribution and a confidence level q
  #' where the t-distribution has df degrees of freedom and 
  #' non-centrality parameter delta.
  #' @param delta Non-centrality parameter
  #' @param test_stat Test statistic at which to evaluate the t-distribution
  #' @param df Degrees of freedom
  #' @param conf_level Confidence level (usually alpha/2 or 1-alpha/2)
  #' @return dif Difference between t-distribution quantile and confidence level
  #' @export

  dif <- pnct(test_stat, df, delta) -  conf_level

  return(dif)
}

getDeltaCI <- function(test_stat, 
                  d,
                  alpha, 
                  n, 
                  interval) {
  #' Confidence intervals for noncentrality parameter of t-distribution
  #'
  #' This function obtains confidence intervals for the non-centrality
  #' parameter of a t-distribution.
  #' @param test_stat Test statistics
  #' @param d Number of parameters in general linear model
  #' @param alpha Significance level
  #' @param n Number of observations in initial study
  #' @param interval Interval within which to search for roots
  #' @return ep Exceedance probability with confidence intervals (vector if cutoff is scalar and matrix otherwise)
  #' @export

  delta_lower <- uniroot(f = tRoot, 
                         interval = interval, 
                         test_stat = test_stat, 
                         df = n - d,
                         conf_level = 1 - alpha / 2)$root

  delta_upper <- uniroot(f = tRoot, 
                         interval = interval, 
                         test_stat = test_stat, 
                         df = n - d,
                         conf_level = alpha / 2)$root

  delta_ci <- c(delta_lower, delta_upper)
  names(delta_ci) <- c("lower", "upper")
  return(delta_ci)
}

getDeltaCI <- Vectorize(getDeltaCI, vectorize.args = "test_stat")


exceedProb <- function(cutoff, 
                  theta_hat, 
                  sd_hat, 
                  d,
                  alpha, 
                  n, 
                  m,
                  interval = c(-100, 100)) {
  #' Confidence intervals for the exceedance probability
  #'
  #' This function obtains confidence intervals for exceedance probability
  #' @param cutoff Cutoff values (scalar or vector)
  #' @param theta_hat Point estimate of parameter of interest
  #' @param sd_hat Estimated standard deviation for parameter of interest
  #' @param d Number of parameters in general linear model
  #' @param alpha Significance level
  #' @param n Number of observations in initial study
  #' @param m Number of observations in replication study
  #' @param interval Interval within which to search for roots
  #' @return ep Exceedance probability with confidence intervals (vector if cutoff is scalar and matrix otherwise)
  #' @export
  #' @examples
  #' library(exceedProb)
  #'
  #' n <- 100
  #' x <- rnorm(n = n)
  #' theta_hat <- mean(x)
  #' sd_hat <- sd(x)
  #' 
  #' cutoff <- seq(from = theta_hat - 0.5, to = theta_hat + 0.5, by = 0.1)
  #' 
  #' exceedProb(cutoff = cutoff, 
  #'            theta_hat = theta_hat, 
  #'            sd_hat = sd_hat, 
  #'            d = 1,
  #'            alpha = 0.05, 
  #'            n = n,
  #'            m = n)

  if (length(theta_hat) != 1) {
    stop("theta_hat must be a scalar")
  }

  if (length(sd_hat) != 1) {
    stop("sd_hat must be a scalar")
  }

  if (length(d) != 1) {
    stop("d must be a scalar")
  }

  if (length(n) != 1) {
    stop("n must be a scalar")
  }

  if (length(m) != 1) {
    stop("m must be a scalar")
  }

  test_stat <- sqrt(n) * (cutoff - theta_hat) / sd_hat

  delta_ci <- getDeltaCI(test_stat = test_stat,
                         d = d,
                         alpha = alpha, 
                         n = n,
                         interval = interval)

  lower <- stats::pnorm(sqrt(m/n) * delta_ci["upper", ], lower.tail = FALSE)
  upper <- stats::pnorm(sqrt(m/n) * delta_ci["lower", ], lower.tail = FALSE)
  point <- stats::pnorm(q = sqrt(m) * (cutoff - theta_hat) / sd_hat, lower.tail = FALSE)

  ep <- data.frame(cutoff = cutoff, point = point, lower = lower, upper = upper)
  return(ep)
}


  
