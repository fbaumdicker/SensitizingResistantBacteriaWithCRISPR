
main <- function(used_cases, observed_copy_numbers, u, dominance_function, number_of_mutations_at_start, type_plasmid_duplication, numerical_number_of_repeats, v, lambda_max, lambda_min, SGV_Type, population_size, competition_ABR_plasmid, competition_CRISPR_plasmid, numerical_approximation_threshold) {
  
  # Packages needed:
  
  library(doParallel)
  library(parallel)
  library(extraDistr)
  library(plyr)
  library(RColorBrewer)
  
  # Directory of file:
  
  main_directory <- dirname(rstudioapi::getSourceEditorContext()$path)
  
  # Source all files (except Standing Genetic Variation, which has to be manually sourced due to Rmd format)
  
  source(paste(main_directory, "/Birth-death-parameters.R", sep = ""))
  source(paste(main_directory, "/First_generation_probabilities.R", sep = ""))
  source(paste(main_directory, "/Numerical_approximation.R", sep = ""))
  source(paste(main_directory, "/Cell_division_plasmid_probabilities.R", sep = ""))
  
  # All observed cases:
  
  cases <- c("Compatible Silencing", "Compatible Cleaving", "Incompatible Cleaving", "Incompatible Silencing", "Compatible No Effect", "Incompatible No Effect")
    
  # Declare Variables:
  
  cell_division_calculated <- list()
  numerical_probabilities <- list()
  rescue_probability_numerical <- list()
  success_rate_simulated <- list()
  distribution_first_generation <- list()
  extinction_probability <- list()
  rescue_probability_SGV <- list()
  rescue_probability_de_novo <- list()
  rescue_probability_total_de_novo <- list()
  
  # Specify all cases:
  
  for (c in used_cases) {
    
    CRISPR_type <- cases[c]
    
    # Calculate the cell division probabilities:
    
    cell_division_calculated_case <- cell_division_plasmid_probabilities(CRISPR_type, observed_copy_numbers, type_plasmid_duplication, dominance_function, competition_ABR_plasmid, competition_CRISPR_plasmid, lambda_max, lambda_min, u)
    
    print(cell_division_calculated_case)
    cell_division_calculated[[length(cell_division_calculated) + 1]] <- cell_division_calculated_case
    
    # Numerical solution for the establishment probabilities:
    
    numerical_probabilities_case <- list()
    
    numerical_probabilities_case <- numerical_establishment_probabilities(observed_copy_numbers, lambda_max, lambda_min, dominance_function, CRISPR_type, cell_division_calculated_case, numerical_approximation_threshold)
    
    numerical_probabilities[[length(numerical_probabilities) + 1]] <- numerical_probabilities_case
    

    # Calculation of the establishment probabilities:
    
    rescue_probability_numerical_case <- rep(0, length(observed_copy_numbers))
    extinction_probability_SGV_case <- rep(0, length(observed_copy_numbers))
    distribution_first_generation_case <- list()
    
    for (j in 1:length(observed_copy_numbers)) {
           
      n <- observed_copy_numbers[j]
      
      ### START FOR CALCULATIONS FOR CELLS WITH FIXED number_of_mutations_at_start pAMRmut AND (n-number_of_mutations) pAMR COPIES (Fig. 3 B)
      
      # Distribution of the first generation after pCRISPR introduction:
      
      first_generation_probabilities_calculated <- first_generation_probabilities(number_of_mutations_at_start, n, CRISPR_type, type_plasmid_duplication, lambda_max, lambda_min, u, dominance_function, competition_CRISPR_plasmid, competition_ABR_plasmid)
      
      # Numerical calculation of the establishment probability:
      
      if (CRISPR_type == "Incompatible Silencing" || CRISPR_type == "Incompatible No Effect") {
        
        for (i in 1:length(first_generation_probabilities_calculated[,1])) {
          if (first_generation_probabilities_calculated[i,1] > 0)
            rescue_probability_numerical_case[j] <- rescue_probability_numerical_case[j] + 
              (1-(1-numerical_probabilities_case[[j]][first_generation_probabilities_calculated[i,1], first_generation_probabilities_calculated[i,2]]) * (1- numerical_probabilities_case[[j]][first_generation_probabilities_calculated[i,3], first_generation_probabilities_calculated[i,4]])) *
              first_generation_probabilities_calculated$probability[i]
        }
        
      } else if (CRISPR_type == "Compatible Silencing") {
        
        if (number_of_mutations_at_start > 0) {
          rescue_probability_numerical_case[j] <- numerical_probabilities_case[[j]][number_of_mutations_at_start]
        } else {
          rescue_probability_numerical_case[j] <- numerical_probabilities_case[[j]][n+1] # = 0
        }
        
      } else if (CRISPR_type == "Compatible Cleaving") {
        
        if (number_of_mutations_at_start > 0) {
          rescue_probability_numerical_case[j] <- numerical_probabilities_case[[j]][n] # = 1 / (1 + lambda_max)
        } else {
          rescue_probability_numerical_case[j] <- 0
        }
        
      } else if (CRISPR_type == "Incompatible Cleaving"){
        
        for (i in 1:length(first_generation_probabilities_calculated[,1])) {
          if (first_generation_probabilities_calculated[i,1] > 0)
            rescue_probability_numerical_case[j] <- rescue_probability_numerical_case[j] + 
              (1-(1-numerical_probabilities_case[[j]][first_generation_probabilities_calculated[i,1]])*(1-numerical_probabilities_case[[j]][first_generation_probabilities_calculated[i,2]])) *
              first_generation_probabilities_calculated$probability[i]
        }
        
      } else if (CRISPR_type == "Compatible No Effect") {
        rescue_probability_numerical_case[j] <- numerical_probabilities_case[[j]][n]
      }
      
      ### START FOR CALCULATIONS FOR CELLS OF A SGV (Fig. 4)
      
      # Calculation of the SGV using Santer's and Uecker's result:
      
      celltypefrequencies_case <- celltypefrequencies_R(v, n, s, SGV_Type)
      
      # Initialize:
      
      if (CRISPR_type == "Incompatible Silencing" || CRISPR_type == "Incompatible No Effect") {
        
        distribution_first_generation_case_copy_number <- matrix(data = 0, nrow = n + 1, ncol = n + 1)
        
      } else {
        
        distribution_first_generation_case_copy_number <- rep(0, n + 1)
        
      }
      
      # The SGV is put together from cells with different number_of_mutations:
      
      for (number_of_mutations in 1:(n+1)) {
        
        first_generation_probabilities_SGV_calculated <- 
          first_generation_probabilities(number_of_mutations, n, CRISPR_type, type_plasmid_duplication, lambda_max, lambda_min, u, dominance_function, competition_CRISPR_plasmid, competition_ABR_plasmid)
        
        if (CRISPR_type == "Incompatible Silencing" || CRISPR_type == "Incompatible No Effect") {
          
          for (i in 1:length(first_generation_probabilities_SGV_calculated[,1])) {
            if (first_generation_probabilities_SGV_calculated[i,1] > 0) {
              distribution_first_generation_case_copy_number[first_generation_probabilities_SGV_calculated[i,1], first_generation_probabilities_SGV_calculated[i,2]] <- distribution_first_generation_case_copy_number[first_generation_probabilities_SGV_calculated[i,1], first_generation_probabilities_SGV_calculated[i,2]] + first_generation_probabilities_SGV_calculated$probability[i] * celltypefrequencies_case[number_of_mutations] * population_size
              distribution_first_generation_case_copy_number[first_generation_probabilities_SGV_calculated[i,3], first_generation_probabilities_SGV_calculated[i,4]] <- distribution_first_generation_case_copy_number[first_generation_probabilities_SGV_calculated[i,3], first_generation_probabilities_SGV_calculated[i,4]] + first_generation_probabilities_SGV_calculated$probability[i] * celltypefrequencies_case[number_of_mutations] * population_size
            }
          }
        } else {
          
          for (i in 1:length(first_generation_probabilities_SGV_calculated[,1])) {
            if (first_generation_probabilities_SGV_calculated[i,1] > 0) {
              distribution_first_generation_case_copy_number[first_generation_probabilities_SGV_calculated[i,1]] <- distribution_first_generation_case_copy_number[first_generation_probabilities_SGV_calculated[i,1]] + first_generation_probabilities_SGV_calculated$probability[i] * celltypefrequencies_case[number_of_mutations] * population_size
              distribution_first_generation_case_copy_number[first_generation_probabilities_SGV_calculated[i,2]] <- distribution_first_generation_case_copy_number[first_generation_probabilities_SGV_calculated[i,2]] + first_generation_probabilities_SGV_calculated$probability[i] * celltypefrequencies_case[number_of_mutations] * population_size
            }
          }
        }
      }
      
      distribution_first_generation_case[[length(distribution_first_generation_case) + 1]] <- distribution_first_generation_case_copy_number
      
      # Calculate the extinction probability of the cells in the SGV:
      
      extinction_probability_SGV_case[j] <- prod((1 - numerical_probabilities_case[[j]]) ^ distribution_first_generation_case_copy_number, na.rm = TRUE)
      
    }
    
    # Calculation of the rescue by de novo mutations only:
  
    if (CRISPR_type == "Compatible Silencing")  {
      rescue_probability_de_novo[[length(rescue_probability_de_novo)+1]] <- rep(0, length(observed_copy_numbers))
      
      for (k in 1:length(observed_copy_numbers)) {
        birth_death_parameter <- birth_death_parameters(observed_copy_numbers[k], lambda_max, CRISPR_type, dominance_function, lambda_min)
        lambda <- birth_death_parameter[[1]]
        mu <- birth_death_parameter[[2]]
        
        poisson_parameter <- v*observed_copy_numbers[k]*population_size*(min(lambda))/(max(mu)-min(lambda))
        
        rescue_probability_de_novo[[length(rescue_probability_de_novo)]][k] <- 1 - sum((1-numerical_probabilities[[c]][[k]][1])^(0:population_size)*dpois(0:population_size,poisson_parameter))
        
       }
    
      } else {rescue_probability_de_novo[[length(rescue_probability_de_novo)+1]] <- rep(0, length(observed_copy_numbers))}
    
    # Save the results:
    
    rescue_probability_numerical[[length(rescue_probability_numerical) + 1]] <- rescue_probability_numerical_case
    distribution_first_generation[[length(distribution_first_generation) + 1]] <- distribution_first_generation_case
    rescue_probability_SGV[[length(rescue_probability_SGV) + 1]] <- 1 - extinction_probability_SGV_case  
    
    }

  
  # Add de novo and SGV probabilities together:
  
  for (c in 1:length(rescue_probability_de_novo)) {
    rescue_probability_total_de_novo[[length(rescue_probability_total_de_novo)+1]] <- 1 - (1 - rescue_probability_de_novo[[c]])*(1-rescue_probability_SGV[[c]])
  }
  
  return(list(rescue_probability_numerical, rescue_probability_SGV, rescue_probability_de_novo, rescue_probability_total_de_novo))
}


