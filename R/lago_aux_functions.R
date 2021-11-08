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
