opt_int <- function(cost, beta, lower, upper, starting.value, pstar, intercept = TRUE, eps = 1.0e-7, max_eval = 3000){
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
  optimum.intervention<-res$solution # Value of the optimal intervention package
  p.obt <- expit(beta, optimum.intervention, intercept) # Value of the obtained outcome goal under the optimal intervention package
  return(list(Optimum_Intervention = optimum.intervention, Obtained_p = p.obt)) # Returning the optimal  intervention package and outcome goal attained under the optimal intervention package
}
