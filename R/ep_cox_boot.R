
CoxEP <- function(data, surv_formula, j, cutoff, n, m, lower_tail) {
  #' Point estimate of the exceedance probability for Cox model parameters
  #' 
  #' This function obtains point estimates for the exceedance probability for Cox model parameters.
  #' @param data survival data (data.frame)
  #' @param surv_formula Survival formula
  #' @param j Index of parameter for which the exceedance probability is obtained
  #' @param cutoff Cutoff values (scalar or vector)
  #' @param n Number of observations in the original study
  #' @param m Number of observations in the replication study (defaults to n if NULL)
  #' @param lower_tail If TRUE, reports lower tail probabilities
  #' @return point (scalar or vector) Point estimate of exceedance probability
  #' @export

  fit0 <- survival::coxph(surv_formula, data = data)
  theta_hat <- stats::coef(fit0)[j]
  sd_hat <- sqrt(n * stats::vcov(fit0)[j, j])

  point <- stats::pnorm(q = sqrt(m) * (cutoff - theta_hat) / sd_hat, lower.tail = lower_tail)

  return(point)
}

exceedProbCoxBoot <- function(data,
                              cox_fit,
                              j, 
                              alpha, 
                              R, 
                              cutoff = NULL, 
                              m = NULL, 
                              lower_tail = FALSE, 
                              sim = "model") {
  #' Bootstrap confidence intervals for the exceedance probability of Cox model parameters
  #' 
  #' This function obtains nonparametric bootstrap percentile confidence for Cox model parameters.
  #' Beta version.
  #' @param data survival data (data.frame)
  #' @param cox_fit (coxph.object) A fitted Cox model
  #' @param j Index of parameter for which the exceedance probability is obtained
  #' @param alpha Significance level
  #' @param R Number of bootstrap resamples
  #' @param cutoff Cutoff values (scalar or vector if supplied, otherwise set to +/- 0.5 of theta_hat)
  #' @param m Number of observations in the replication study (defaults to n if NULL)
  #' @param lower_tail If TRUE, reports lower tail probabilities; otherwise reports upper tail probabilities
  #' @param sim type of simulation, input to boot::censboot
  #' @return ep Exceedance probability with confidence intervals
  #' @export
  #' @examples
  #' library(exceedProb)
  #' library(survival)
  #' 
  #' # Cox model -------------------------------------------------------
  #'
  #' # Simulate exponential data
  #' n <- 50
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

  #' time <- pmin(event_time, censor_time)
  #' event <- time == event_time
  #' surv_data = data.frame(time = time, event = event, group = tx_indicator)
  #'
  #' # Fit Cox model and get bootstrap percentile confidence intervals for the exceedance probability 
  #' # with model-based resampling (see documentation for boot::censboot)
  #' cox_fit <- coxph(Surv(time, event) ~ group, data = surv_data)
  #' ep <- exceedProbCoxBoot(data = surv_data,
  #'                         cox_fit = cox_fit,
  #'                         j = 1, 
  #'                         alpha = 0.05, 
  #'                         R = 500)
  #'
  #' # Plot results
  #' with(ep, plot(cutoff, point, type = "l"))
  #' with(ep, lines(cutoff, lower, lty = 2))
  #' with(ep, lines(cutoff, upper, lty = 2))

  event_chr <- as.character(cox_fit$formula[[2]])[3]

  # Todo(bsegal): revise code below so this warning is not necessary
  if(!is.logical(eval(parse(text = event_chr), data))) {
    stop(paste0("The event indicator must evaluate to a logical; ", event_chr, " does not."))
  }

  # TODO(bsegal): Add other checks to ensure input is valid

  km_surv <- survival::survfit(cox_fit)

  # pull out elements of model formula
  time_chr <- as.character(cox_fit$formula[[2]])[2]
  covariate_chr <- as.character(cox_fit$formula[[3]])

  # create formula for reverse KM
  time_revised_chr <- paste0(time_chr, " - 0.001*(", event_chr, ")")
  event_reverse_chr <- paste0("!(", event_chr, ")")
  rev_km_formula <- stats::as.formula(paste0("Surv(", time_revised_chr, ", ", event_reverse_chr, ") ~ ", covariate_chr))

  # assuming censoring distribution explained by same covariates as survival distribution
  km_cens <- survival::survfit(rev_km_formula, data = data)

  # using number of rows to scale variance,
  # assuming number of events is proportional to number of rows
  n <- nrow(data)

  if (is.null(m)) {
    m <- n
  }

  theta_hat <- stats::coef(cox_fit)[j]
  sd_hat <- sqrt(n * stats::vcov(cox_fit)[j, j])

  if(is.null(cutoff)) {
    cutoff <- seq(from = theta_hat - 1, to = theta_hat + 1, by = 0.1)
  }

  point <- stats::pnorm(q = sqrt(m) * (cutoff - theta_hat) / sd_hat, lower.tail = lower_tail)

  boot_out <- boot::censboot(data = data,
                             statistic = CoxEP,
                             R = R,
                             F.surv = km_surv,
                             G.surv = km_cens,
                             sim = sim,
                             cox = cox_fit,
                             surv_formula = cox_fit$formula,
                             j = j,
                             cutoff = cutoff,
                             n = n,
                             m = m,
                             lower_tail = lower_tail)

  ci <- apply(boot_out$t, 2, stats::quantile, probs = c(1 - alpha / 2, alpha / 2), type = 6)
  ep <- data.frame(cutoff = cutoff, point = point, lower = ci[2, ], upper = ci[1, ])

  return(ep)
}
