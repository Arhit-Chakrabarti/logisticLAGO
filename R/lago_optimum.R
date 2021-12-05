#' Optimum intervention under the LAGO design
#' @description The first goal of a Learn-As-You-Go (LAGO) study is to identify the optimal intervention package while minimizing the cost of the intervention package subject to the probability of a desired binary outcome being above a given threshold. The current LAGO design considers a logistic regression model with linear cost coefficients.
#'  This function takes as input the vector of linear cost coefficients, the true/estimated beta values, the lower and upper limits for the components of the intervention package, the desired outcome goal and a starting value of the intervention package for the optimization algorithm and calculates the optimal intervention package.
#'
#' @param cost a vector of linear cost coefficients
#' @param beta a vector of true/estimated beta values
#' @param lower a vector providing the values of lower limits for the components of the intervention package
#' @param upper a vector providing the values of upper limits for the components of the intervention package
#' @param starting.value a vector of starting values for the components of the intervention package
#' @param pstar the desired outcome goal
#' @param intercept a logical argument to include intercept in the model or not (default = TRUE)
#' @param eps the desired level of tolerance (default = 1.0e-7)
#' @param max_eval the maximum number of iterations to perform (default = 3000)
#'
#' @return The returned value is a named list containing
#' \item{Optimum_Intervention}{the optimum intervention package}
#' \item{Obtained_p}{the outcome goal under the optimal intervention package}
#' @export
#'
#' @examples
#' # Defining vector of starting values for the algorithm
#' x.init = c(2.5, 12.5, 7)
#' # Defining vector of lower limits for the components
#' x.l = c(1, 10, 2)
#' # Defining vector of upper limits for the components
#' x.u = c(4, 15, 15)
#' # Defining vector of linear cost coefficients
#' cost_lin = c(1, 8, 2.5)
#' p_bar = 0.9 # Defining the desired outcome goal
#' # True beta values for the study
#' beta = c(log(0.05), log(1.2), log(1.1), log(1.3))
#' ## Running the LAGO optimization algorithm
#' opt_lago = opt_int(cost = cost_lin, beta = beta, lower = x.l,
#'                    upper = x.u, pstar = p_bar, starting.value = x.init)
opt_int <- function(cost, beta, lower, upper, starting.value, pstar, intercept = TRUE, eps = 1.0e-7, max_eval = 3000){
  # Checks
  if(length(upper) != length(lower)) stop("lower and upper limits of the components are not of the same length") # Length of upper should be equal to length of lower ranges of the components
  if(length(upper) != length(starting.value)) stop("length of starting.value and upper and lower do not match") # Provided ranges of the components should have the same dimension as that of the intervention package itself
  if (any(lower < 0)) stop("The intervention must have non-negative values only") # Non negative interventions only allowed
  if (any(lower >= upper)) stop("Upper limits of the intervention package must be larger than corresponding lower limits")  # Components cannot have lower ranges more than or equal to the upper values
  if(intercept == TRUE & length(beta) != (length(starting.value) + 1)){
    stop("Please provide the correspodning beta0 value if the model should include the intercept. Else please change intercept to FALSE")
  }
  if(intercept == FALSE & length(beta) != length(starting.value)){
    stop("Please check the dimension of provided beta value")
  }
  if(length(cost) != length(starting.value)) stop("length of cost vector and the length of intervention package do not match")  # Cost vector should be equal to the number of components in the package
  if(pstar < 0 | pstar > 1) stop("desired probability of success can be between 0 and 1") # Probability should be between 0 and 1
  #Defining Objective Function
  eval_f <- function(x) {
    obj<- sum(cost * x) # Total cost of the intervention package
    grad <- cost # Calculating the gradient of the objective function
    return(list("objective" = obj, "gradient" = grad)) # Returning the list of value of objective function and gradient
  }
  # Take out the intercept from the model for calculating the Jacobian of the inequality
  if(intercept == TRUE){
    beta_without_intercept <- beta[-1] # Beta value without the intercept
  }else{
    beta_without_intercept <- beta
  }
  #Defining the constraint
  eval_g_ineq <- function(x) {
    prob = expit(beta, x, intercept) # Calculating the probability of success
    constr <- (pstar - prob) # Defining the outcome goal constraint
    grad = - ((prob^2) * exp_my(beta = - beta, x, intercept) * beta_without_intercept) # Calculating gradient of inequality consraint
    return(list( "constraints" = constr, "jacobian" = grad )) # Returning the list of value of constraint function and gradient
  }

  ## Defining the Algorithm Options for optimizing using non-linear optimization techniques
  local_opts <- list( "algorithm" = "NLOPT_LD_MMA", # Using Method of Moving Asymptotes
                      "xtol_rel" = eps )
  opts <- list( "algorithm" = "NLOPT_LD_AUGLAG_EQ", # Augmented Lagrangian method ( Or "NLOPT_LD_AUGLAG")
                "xtol_rel" = eps,
                "maxeval" = max_eval,
                "local_opts" = local_opts )
  x0 <- starting.value    # Starting values of Intervention for the optimization algorithm
  lb <- lower	 # Lower Bounds
  ub <- upper	 # Upper Bounds

  # Optimization using non-linear optimization
  res <- nloptr::nloptr(x0 = x0, eval_f = eval_f, lb = lb, ub = ub, eval_g_ineq = eval_g_ineq, opts = opts)
  optimum.intervention <- res$solution # Value of the optimal intervention package
  p.obt <- expit(beta, optimum.intervention, intercept) # Value of the obtained outcome goal under the optimal intervention package

  return(list(Optimum_Intervention = optimum.intervention, Obtained_p = p.obt)) # Returning the optimal  intervention package and outcome goal attained under the optimal intervention package
}
