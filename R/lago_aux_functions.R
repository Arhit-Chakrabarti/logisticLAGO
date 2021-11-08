## Defining the Indicator function
Indicator <- function(condition){
  ifelse(condition, 1, 0)
}
## Defining the expit function i.e. 1/(1 + exp(-\beta^T X)) which takes as input, beta, x and a logical input to include intercept or not
expit <- function(beta, x, intercept = TRUE){
  # If intercept is TRUE i.e. in supplied parameter beta an intercept term is provided
  if(intercept == TRUE){
    u <- c(1, x) # Create a new vector with component 1 and the other components as the x values
  }else{ # If in supplied beta the component for intercept is not supplied
    u <- x # Create a new vector same as the x values
  }
  linear_comp <- sum(beta * u) # Calculate the value of the linear predictor \beta^T X
  return(1 / (1 + exp(- linear_comp))) # Return the value of the expit function
}
## Defining the expit function which takes as input the value of the linear predictor as input
expit_linpred <-function(x){
  return(1/(1+exp(-x)))
}
## Defining the my exponential function i.e. exp(-\beta^T X)) which takes as input, beta, x and a logical input to include intercept or not
exp_my <- function(beta, x, intercept = TRUE){
  # If intercept is TRUE i.e. in supplied parameter beta an intercept term is provided
  if(intercept == TRUE){
    u <- c(1, x) # Create a new vector with component 1 and the other components as the x values
  }else{# If in supplied beta the component for intercept is not supplied
    u <- x # Create a new vector same as the x values
  }
  linear_comp <- sum(beta * u) # Calculate the value of the linear predictor \beta^T X
  return(exp(linear_comp)) # Return the value of my exp function
}
## Defining the function to jitter observations for the simulation codes
jitter_my <- function(x, n, variation){# Takes as input the vector of values x, sample size n and the variation parameter which controls the amount of jittering induced by the function about the initial x values
  x_jitter = x + (variation * x * t(matrix(stats::rnorm(n * length(x), mean = 0, sd = 1), nrow = n)))
  return(t(x_jitter))
}
## Defining the truncate function to truncate the actual interventions after jittering the starting/optimum intervention at any stage
truncate <- function(x, lower.limit, upper.limit){
  # If any value x is below the lower range for the corresponding component, replace the value with the lower limit
  first_part = lower.limit * Indicator(t(x) < lower.limit)
  # If the value of x lies between the lower and upper limits then return the corresponding value
  second_part = t(x) * Indicator(lower.limit < t(x)) * Indicator(t(x) < upper.limit)
  # If any value x is above the upper range for the corresponding component, replace the value with the upper limit
  third_part = upper.limit * Indicator(t(x) > upper.limit)
  truncated = first_part + second_part + third_part # Adding the three cases together
  return(t(truncated)) # Returning the vector of truncated values
}


