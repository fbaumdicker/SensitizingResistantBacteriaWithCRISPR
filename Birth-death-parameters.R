
# Calculation of the birth and death rates for a cell with given copy number and CRISPR type.

# OUTCOME: The birth and death rates.


birth_death_parameters <- function(n, lambda_max, CRISPR_type, dominance_function, lambda_min) {
  
  if (CRISPR_type == "Incompatible No Effect") { # Nothing is silenced! 
      
    lambda <- matrix(data = 1 + lambda_max, nrow = n+1, ncol = n+1)
    mu <- matrix(data = 1, nrow = n+1, ncol = n+1) 
    
    if (dominance_function == "dominant") {
      
      lambda[n, n+1] <- 1 - lambda_min
      
    } else if (dominance_function == "recessive") {
      
      lambda <- matrix(data = 1 - lambda_min, nrow = n+1, ncol = n+1)
      
      lambda[n+1,] <- 1 + lambda_max
      
    } else if (dominance_function == "linear") {
      
      for (col in 1:(n+1)) {
        lambda[, col] <- seq(from = 1 + lambda_max, to = 1 - lambda_min, length.out = n + 1) 
      }
      
      lambda <- rbind(lambda[2:(n+1),], lambda[1,])
      
    } 
      
  } else if (CRISPR_type == "Incompatible Silencing") {
    
    lambda <- matrix(data = 1 + lambda_max, nrow = n+1, ncol = n+1)
    mu <- matrix(data = 1, nrow = n+1, ncol = n+1) 
    
    if (dominance_function == "dominant") {
      
      lambda[1:n,n+1] <- 1 - lambda_min
      
    } else if (dominance_function == "recessive") {
      
      lambda[1:n,] <- 1 - lambda_min
      
    } else if (dominance_function == "linear") {
      
      for (row in 1:n) {
        lambda[row, ] <- seq(from = 1 - lambda_min, to = 1 + lambda_max, length.out = n + 1)
      }
      lambda <- cbind(lambda[,2:(n+1)], lambda[,1])
      
    }
  } else if (CRISPR_type == "Compatible Cleaving" || CRISPR_type == "Incompatible Cleaving" || CRISPR_type == "Compatible Silencing") {
      
      lambda <- rep(1 + lambda_max, n + 1)
      mu <- rep(1, n+1)  
      
      if(dominance_function == "dominant") {
        
        lambda[n+1] <- 1 - lambda_min
        
      } else if (dominance_function == "recessive") {
        
        lambda[n+1] <- 1 - lambda_min
        if (n > 1) lambda[1:(n-1)] <- 1 - lambda_min
        
      } else if (dominance_function == "linear") {
        
        lambda <- seq(from = 1 - lambda_min, to = 1 + lambda_max, length.out = n + 1)
        lambda <- c(lambda[2:(n+1)], lambda[1])
        
      } 
      
  } else if (CRISPR_type == "Compatible No Effect") {
    
    lambda <- rep(1 + lambda_max, n + 1)
    mu <- rep(1, n+1) 
    
  } else {return("Error (wrong type)")}
  
  return(list(lambda, mu))
  
}

