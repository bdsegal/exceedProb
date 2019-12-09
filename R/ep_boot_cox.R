
fitCox <- function(data, surv_formula, j, cutoff, n, m, lower_tail) {
  #' Point estimate of the exceedance probability for Cox model parameters
  #' 
  #' This function obtains point estimates for the exceedance probability for Cox model parameters.
  #' @param data survival data (data.frame)
  #' @param surv_formula Survival formula (character), e.g. "Surv(time, event) ~ group"
  #' @param j  Index of parameter for which the exceedance probability is obtained
  #' @param cutoff Cutoff values (scalar or vector)
  #' @param n Number of observations in the original study
  #' @param m Number of observations in the replication study (defaults to n if NULL)
  #' @param lower_tail If TRUE, reports lower tail probabilities
  #' @return point (scalar or vector) Point estimate of exceedance probability
  #' @export

  fit0 <- survival::coxph(stats::as.formula(paste0("survival::", surv_formula)), data = data)
  theta_hat <- stats::coef(fit0)[j]
  sd_hat <- sqrt(n * stats::vcov(fit0)[j, j])

  point <- stats::pnorm(q = sqrt(m) * (cutoff - theta_hat) / sd_hat, lower.tail = lower_tail)

  return(point)
}

exceedProbCoxBoot <- function(data, surv_formula, j, cutoff = NULL, m = NULL, alpha = 0.05, B = 1000, lower_tail = FALSE) {
  #' Bootstrap confidence intervals for the exceedance probability of Cox model parameters
  #' 
  #' This function obtains nonparametric bootstrap percentile confidence for Cox model parameters.
  #' Beta version.
  #' @param data survival data (data.frame)
  #' @param surv_formula Survival formula (character), e.g. "Surv(time, event) ~ group"
  #' @param j Index of parameter for which the exceedance probability is obtained
  #' @param cutoff Cutoff values (scalar or vector if supplied, otherwise set to +/- 0.5 of theta_hat)
  #' @param alpha Significance level
  #' @param m Number of observations in the replication study (defaults to n if NULL)
  #' @param B Number of bootstrap resamples
  #' @param lower_tail If TRUE, reports lower tail probabilities
  #' @return ep Exceedance probability with confidence intervals
  #' @export
  #' @examples
  #' library(exceedProb)
  #' 
  #' # Cox model -------------------------------------------------------
  #'
  #' # Simulate exponential data
  #' n <- 100
  #' m <- 100
  #' baseline_hazard <- 1
  #' theta <- 0.4
  #' p_censor <- 0.3
  #' prop_tx <- 0.5
  #'
  #' tx_indicator = rbinom(n = n, size = 1, prob = prop_tx)
  #' event_rate <- baseline_hazard * exp(theta * tx_indicator)
  #' censor_rate <- event_rate * p_censor / (1 - p_censor)
  #'
  #' event_time <- rexp(n = n, rate = event_rate)
  #' censor_time <- rexp(n = n, rate = censor_rate)
  #'
  #' time <- pmin(event_time, censor_time)
  #' event <- time == event_time
  #'
  #' surv_data = data.frame(time = time, event = event, group = tx_indicator)
  #'
  #' # Get exceedance probability with nonparametric bootstrap percentile confidence intervals
  #' ep <- exceedProbCoxBoot(data = surv_data,
  #'                         surv_formula = "Surv(time, event) ~ group",
  #'                         j = 1,
  #'                         B = 1000)

  # TODO(bsegal): Add checks to ensure input is valid

  n <- nrow(data)

  if(is.null(m)) {
    m <- n
  }

  fit0 <- survival::coxph(stats::as.formula(paste0("survival::", surv_formula)), data = data)
  theta_hat <- stats::coef(fit0)[j]
  sd_hat <- sqrt(n * stats::vcov(fit0)[j, j])

  if(is.null(cutoff)) {
    cutoff <- seq(from = theta_hat - 0.5, to = theta_hat + 0.5, by = 0.1)
  }

  point <- stats::pnorm(q = sqrt(m) * (cutoff - theta_hat) / sd_hat, lower.tail = lower_tail)

  boot_out <- matrix(NA, ncol = length(cutoff), nrow = B)

  for (b in 1:B) {
    data_star <- data[sample(n, size = n, replace = TRUE), ]
    boot_out[b, ] <- fitCox(data = data_star,
                            surv_formula = surv_formula,
                            j = j,
                            cutoff = cutoff,
                            n = n,
                            m = m,
                            lower_tail = lower_tail)
  }

  ci <- apply(boot_out, 2, stats::quantile, probs = c(alpha / 2, 1 - alpha / 2))

  ep <- data.frame(cutoff = cutoff, point = point, lower = ci[1, ], upper = ci[2, ])
  rownames(ep) <- NULL

  return(ep)
}
