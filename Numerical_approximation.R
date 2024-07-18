
# Computation of the establishment probabilities / the probability to go extinct:

# OUTCOME: Extinction probabilities for each cell type


# 1. Definition of the probability generation function. 

# Variables:
### Q_input: The variable (input) of the probability generating function.
### cell_division_calculated: Probabilities for each combination of daughter cells
### lambda: Birth rates
### mu: Death rates
### n: copy number
### k: Position of the copy number n in the associated cell_division_calculated list.

probability_generating_function <- function(Q_input, CRISPR_type, cell_division_calculated, lambda, mu, n, k) {
  
  if (CRISPR_type %in% c("Incompatible Cleaving", "Compatible Cleaving", "Compatible Silencing", "Compatible No Effect")) {
    Q_output <- c(rep(0,n),0)
    
    for (ABR_plasmids_in_cell in 1:(n+1)) {
      
      Q_output[ABR_plasmids_in_cell] <- Q_output[ABR_plasmids_in_cell] + sum(Q_input[cell_division_calculated[[k]][[ABR_plasmids_in_cell]][,1]]*Q_input[cell_division_calculated[[k]][[ABR_plasmids_in_cell]][,2]]*cell_division_calculated[[k]][[ABR_plasmids_in_cell]][,3])
      
      Q_output[ABR_plasmids_in_cell] <- 1/(lambda[ABR_plasmids_in_cell]+mu[ABR_plasmids_in_cell])* (mu[ABR_plasmids_in_cell] + lambda[ABR_plasmids_in_cell]*Q_output[ABR_plasmids_in_cell])
      # Notion: u is not needed in the formula, since we only determine if the cell dies or it divides itself. u is then included only in the calculations only if the cell divides.
    }
  }
  
  else if (CRISPR_type %in% c("Incompatible Silencing", "Incompatible No Effect")) {
    Q_output <- matrix(data = 0, nrow = n+1, ncol = n+1)
    
    for (CRISPR_plasmids_in_cell in 1:(n+1)) {
      for (mutant_ABR_plasmids_in_cell in 1:(n+1)) {
        
        if (CRISPR_plasmids_in_cell%%(n+1) + mutant_ABR_plasmids_in_cell%%(n+1) <= n) {
          Q_output[CRISPR_plasmids_in_cell, mutant_ABR_plasmids_in_cell] <- Q_output[CRISPR_plasmids_in_cell, mutant_ABR_plasmids_in_cell] + sum(Q_input[as.matrix(cell_division_calculated[[k]][[CRISPR_plasmids_in_cell]][[mutant_ABR_plasmids_in_cell]][,1:2])]*Q_input[as.matrix(cell_division_calculated[[k]][[CRISPR_plasmids_in_cell]][[mutant_ABR_plasmids_in_cell]][,3:4])]*cell_division_calculated[[k]][[CRISPR_plasmids_in_cell]][[mutant_ABR_plasmids_in_cell]][,5])
          
          Q_output[CRISPR_plasmids_in_cell, mutant_ABR_plasmids_in_cell] <- 1/(lambda[CRISPR_plasmids_in_cell, mutant_ABR_plasmids_in_cell] + mu[CRISPR_plasmids_in_cell, mutant_ABR_plasmids_in_cell])* (mu[CRISPR_plasmids_in_cell, mutant_ABR_plasmids_in_cell] + lambda[CRISPR_plasmids_in_cell, mutant_ABR_plasmids_in_cell]*Q_output[CRISPR_plasmids_in_cell, mutant_ABR_plasmids_in_cell])
          # Notion: u is not needed in the formula, since we only determine if the cell dies or it divides itself. u is then included only in the calculations only if the cell divides.
        } else Q_output[CRISPR_plasmids_in_cell, mutant_ABR_plasmids_in_cell] <- NaN
      }
    }
  }
  
  else {print("Wrong CRISPR type")}
  
  
  return(Q_output)
}


# 2. Numerical calculation of the establishment probabilities. 

# Variables:
### observed_copy_numbers: The observed copy numbers n.
### s_max: The maximum fitness.
### s_0: The minimum fitness.
### dominance_function: The fitness function (dominant, recessive or linear).
### CRISPR_type: The CRISPR type and the incompatability group.
### cell_division_calculated: Probabilities for each combination of daughter cells.
### numerical_approximation_threshold: The numerical approximation runs until the change in value is less than this parameter.

numerical_establishment_probabilities <- function(observed_copy_numbers, s_max, s_0, dominance_function, CRISPR_type, cell_division_calculated, numerical_approximation_threshold) {
  
  computation_result <- list()
  
  for (j in 1:length(observed_copy_numbers)) {
    
    # Get copy number and birth and death rates:
    
    n <- observed_copy_numbers[j]
    
    birth_death_parameter <- birth_death_parameters(n, s_max, CRISPR_type, dominance_function, s_0)
    lambda <- birth_death_parameter[[1]]
    mu <- birth_death_parameter[[2]]
    
    # Numerical approximation:

    if (CRISPR_type %in% c("Incompatible Silencing", "Incompatible No Effect")) {
      
      # Initial value 0:
      
      z <- matrix(data = 0, nrow = n+1, ncol = n+1)
      
      # First iteration:
      
      z_new <- probability_generating_function(z, CRISPR_type, cell_division_calculated, lambda, mu, n, j)
      
      # While the change in value is bigger than numerical_approximation_threshold, the approximation is continued.
      
      while (max(abs(z_new - z), na.rm = TRUE) >= numerical_approximation_threshold) {
        z <- z_new
        z_new <- probability_generating_function(z, CRISPR_type, cell_division_calculated, lambda, mu, n, j)
      }
      
      computation_result[[length(computation_result) + 1]] <- 1-z
      
    } else if (CRISPR_type %in% c("Incompatible Cleaving", "Compatible Cleaving", "Compatible Silencing", "Compatible No Effect")) {
      
      # Initial value 0:
      
      z <- rep(0, n + 1)
      
      # First iteration:
      
      z_new <- probability_generating_function(z, CRISPR_type, cell_division_calculated, lambda, mu, n, j)
      
      # While the change in value is bigger than numerical_approximation_threshold, the approximation is continued.
      
      while (max(abs(z_new - z)) >= numerical_approximation_threshold) {
        z <- z_new
        z_new <- probability_generating_function(z, CRISPR_type, cell_division_calculated, lambda, mu, n, j)
      }
      
      computation_result[[length(computation_result) + 1]] <- 1-z
    
      }
  }
  
  return(computation_result)
  
}


