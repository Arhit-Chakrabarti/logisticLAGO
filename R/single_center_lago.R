#' Optimum package under the single center LAGO design
#'
#' @param x0 a vector of starting values for the components of the intervention package
#' @param lower a vector providing the values of lower limits for the components of the intervention package
#' @param upper a vector providing the values of upper limits for the components of the intervention package
#' @param nstages the number of stages of the LAGO design K (>=2)
#' @param beta.true a vector of true i.e. population beta values
#' @param sample.size sample size (n) in each of the K stages
#' @param icc the expected variation among subjects within the center at any stage
#' @param cost.vec a vector of per unit linear cost coefficients
#' @param prob the desired outcome goal
#' @param intercept a logical argument to include intercept in the model or not (default = TRUE)
#' @param B the number of simulations (default = 100)
#'
#' @return
#' @export
#'
#' @examples
sc_lago <- function(x0, lower, upper, nstages, beta.true, sample.size, icc, cost.vec, prob, intercept = TRUE, B = 100){
  # Initializing objects for use in the simulation
  response_clubbed <- NULL # Response
  actual_intervention_clubbed <- NULL # Actual intervention
  xopt <- array(NA, dim = c(nstages, length(x0), B)) # Array to store the optimal intervention
  p.opt.hat <- matrix(NA, nrow = nstages, ncol = B) # Matrix to store the obtained outcome goal
  power <- rep(NA, times = B) # To store the power
  for(b in 1:B){
    x.start = x0 # Initial value of x at each simulation
    for(k in 1:nstages){
      actual_intervention = jitter_my(x.start, n = sample.size, variation = icc)
      actual_intervention = truncate(x = actual_intervention, lower.limit = lower, upper.limit = upper)
      actual_intervention_clubbed = rbind(actual_intervention_clubbed, actual_intervention)

      actual_intervention_full <- cbind(1, actual_intervention)
      lin_pred = t(beta.true * t(actual_intervention_full))
      response <- rep(NA, sample.size)
      for(i in 1:sample.size){
        response[i] <- stats::rbinom(n = 1, size = 1, prob = expit_linpred(x = sum(lin_pred[i, ])))
      }
      response_clubbed <- c(response_clubbed, response)
      fit <- stats::glm(response_clubbed ~ actual_intervention_clubbed, family = "binomial")
      beta_hat = unname(fit$coefficients)
      opt_lago = logisticLAGO::opt_int(cost = cost.vec, beta = beta_hat, lower = lower, upper = upper, pstar = prob, starting.value = x.start, intercept)
      xopt[k, ,b] = opt_lago$Optimum_Intervention
      x.start = xopt[k, ,b]
      p.opt.hat[k, b] = opt_lago$Obtained_p
    }
    # Power calculation
    if(intercept == TRUE){
      test_stat <- sum(beta_hat[ - 1] * solve(stats::vcov(fit)[ - 1, - 1], beta_hat[ - 1]))
    }else{
      test_stat <- sum(beta_hat * solve(stats::vcov(fit), beta_hat))
    }
    p.value <- stats::pchisq(test_stat, df = (length(beta_hat) - 1), lower.tail = FALSE) # Calculate p-value
    power[b] <- ifelse(p.value <= 0.05, 1, 0) # Calculate power
  }
  opt.x <- as.data.frame(apply(xopt, MARGIN = c(1, 2), FUN = mean))
  opt.p <- apply(p.opt.hat, MARGIN = 1, FUN = mean)
  power.obt <- mean(power)
  name.col = 0
  for(i in 1:length(x0)){name.col[i] <- (paste0("x",i))}
  name.row = 0
  for(i in 1:nstages){name.row[i] <- (paste0("Stage",i))}
  colnames(opt.x) <- name.col; rownames(opt.x) <- name.row; names(opt.p) <- name.row

  return(list(xopt = opt.x, p.opt.hat = opt.p, power = power.obt))
}
