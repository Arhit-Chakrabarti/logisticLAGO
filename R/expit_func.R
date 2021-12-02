#' Success probabilities under the logistic regression model
#' @description This function calculates the success probabilities for the logistic regression model for a given vector of covariates and a vector of parameter values. The success probability for a given vector of covariates x and parameter vector beta is calculated as 1/(1 + exp(-beta' X))
#' @param beta a vector of true/estimated values of the parameters for a logistic regression model. If intercept = TRUE is specified, the first component should include the value of the intercept term
#' @param x a vector of covariates under the logistic regression model
#' @param intercept a logical argument to specify if the intercept term is included in the model or not. If TRUE, then beta should include the intercept term (default = TRUE)
#'
#' @return This function returns the calculated probability of success under the logistic regression model for the given covariate vector and vector of parameters.
#' @export
#'
#' @examples
#' x = c(2.5, 12.5, 7)
#' beta.true  = c(log(0.05), log(1.2), log(1.1), log(1.3))
#' # Default intercept = TRUE
#' prob <- expit(beta.true, x) # Calculating the probabilities
#'
#' beta.true  = c(log(1.2), log(1.1), log(1.3))
#' # intercept = FALSE
#' prob <- expit(beta.true, x, intercept = FALSE) # Calculating the probabilities
expit <- function(beta, x, intercept = TRUE){
  if(intercept == TRUE & length(beta) != (length(x) + 1)){
    stop("argument specifies intercept == TRUE, but beta does not include intercept term")
  }
  if(intercept == FALSE & length(beta) != length(x)){
    stop("The length of x and beta do not match")
  }
  # If intercept is TRUE i.e. in supplied parameter beta an intercept term is provided
  if(intercept == TRUE){
    u <- c(1, x) # Create a new vector with component 1 and the other components as the x values
  }else{ # If in supplied beta the component for intercept is not supplied
    u <- x # Create a new vector same as the x values
  }
  linear_comp <- sum(beta * u) # Calculate the value of the linear predictor \beta^T X
  return(1 / (1 + exp(- linear_comp))) # Return the value of the expit function
}
