

# Calculation of all possible cell type combinations of daughter cells, given number of plasmids in the mother cell and the type of duplication.
# The probability of dying instead of replicating is included.

# OUTPUT: Probability for all possible daughter cell type combinations.

cell_division_plasmid_probabilities <- function(CRISPR_type, observed_copy_numbers, type_plasmid_duplication, dominance_function, competition_AMR_plasmid, competition_CRISPR_plasmid, lambda_max, lambda_min, u) {
  
  if (CRISPR_type == "Incompatible Silencing" || CRISPR_type == "Incompatible No Effect") {
    
    # Create a list for the result. The first entries of the list is the copy number, the second is the number of CRISPR plasmids and the third is the number of mutated AMR plasmids. The third is the number of additional CRISPR plasmids, gained additionally between the cell division.
    
    result <- vector("list", length = length(observed_copy_numbers))
    
    for (l in 1:length(observed_copy_numbers)) {
      n <- observed_copy_numbers[l]
      result[[l]] <- vector("list", length = n + 1)
      
      for (i in 1:(n + 1)) {
        result[[l]][[i]] <- vector("list", length = n + 1)
      }
    }
    
    # Calculate the probabilities for the daughter cells for all mother cells and observed copy numbers:
    
    for (i in 1:length(observed_copy_numbers)) {
      
      n <- observed_copy_numbers[i] # Copy number
      
      birth_death_parameter <- birth_death_parameters(n, lambda_max, CRISPR_type, dominance_function, lambda_min) # Birth and death parameters for the specific cell type
        
      for (mother_cell_CRISPR_plasmids in 1:(n+1)) {
        for (mother_cell_mutated_AMR_plasmids in c(1:(n-mother_cell_CRISPR_plasmids%%(n+1)), n+1)) {
          
          if (mother_cell_CRISPR_plasmids == n) mother_cell_mutated_AMR_plasmids = n+1
          
          # Probability to gain more CRISPR encoding plasmids:
          
          probability_u <- u / (birth_death_parameter[[1]][mother_cell_CRISPR_plasmids, mother_cell_mutated_AMR_plasmids] + u)
          
          # Creation of a matrix with all possible combinations of daughter cells:
          
          probabilities_plasmids_in_daugther_cells <- expand.grid(
            CRISPR_first_daugther = 1:(n+1),
            mutant_AMR_first_daugther = 1:(n+1),
            CRISPR_second_daugther = 1:(n+1),
            mutant_AMR_second_daugther = 1:(n+1)
          )
          
          # Calculation of the probabilities of the amount of plasmids after plasmid duplication (before cell division): 
          
          probabilities_plasmids_in_daugther_cells <- subset(probabilities_plasmids_in_daugther_cells, probabilities_plasmids_in_daugther_cells$CRISPR_first_daugther%%(n+1) + probabilities_plasmids_in_daugther_cells$mutant_AMR_first_daugther%%(n+1)  <= n &
                                                               probabilities_plasmids_in_daugther_cells$CRISPR_second_daugther%%(n+1) + probabilities_plasmids_in_daugther_cells$mutant_AMR_second_daugther%%(n+1)  <= n)
          
          probabilities_plasmids_in_daugther_cells$CRISPR_parent_after_doubling <- probabilities_plasmids_in_daugther_cells$CRISPR_first_daugther%%(n+1) + probabilities_plasmids_in_daugther_cells$CRISPR_second_daugther%%(n+1)
          probabilities_plasmids_in_daugther_cells$mutant_AMR_parent_after_doubling <- probabilities_plasmids_in_daugther_cells$mutant_AMR_first_daugther%%(n+1) + probabilities_plasmids_in_daugther_cells$mutant_AMR_second_daugther%%(n+1)
          
          probabilities_plasmids_in_daugther_cells$CRISPR_parent_after_doubling <- probabilities_plasmids_in_daugther_cells$CRISPR_parent_after_doubling + (2*n+1)*(probabilities_plasmids_in_daugther_cells$CRISPR_parent_after_doubling == 0)
          probabilities_plasmids_in_daugther_cells$mutant_AMR_parent_after_doubling <- probabilities_plasmids_in_daugther_cells$mutant_AMR_parent_after_doubling + (2*n+1)*(probabilities_plasmids_in_daugther_cells$mutant_AMR_parent_after_doubling == 0)
          
          probabilities_plasmids_in_daugther_cells$probability <- 0
          
          for (additional_CRISPR_plasmids in 0:n) {
  
            probabilities_plasmids_after_duplication <- matrix(data = 0, nrow = 2*n+1, ncol = 2*n+1)
            
            plasmids_in_cell <- c(mother_cell_CRISPR_plasmids%%(n+1) + additional_CRISPR_plasmids, mother_cell_mutated_AMR_plasmids%%(n+1), 
                                  n - mother_cell_CRISPR_plasmids%%(n+1) - mother_cell_mutated_AMR_plasmids%%(n+1))
            
            
            if (type_plasmid_duplication == "random") {
              
              # Polya's urn model: Only one plasmid is replicated at a time. The plasmid is chosen randomly. The duplicated plasmid is added to the pool.
              
              for (new_CRISPR_plasmids in 0:(n-additional_CRISPR_plasmids)) {
                for (new_mutant_AMR_plasmids in 0:(n-additional_CRISPR_plasmids)) {
                  if(new_CRISPR_plasmids + plasmids_in_cell[1] > 0) row <- new_CRISPR_plasmids + plasmids_in_cell[1] else row <- 2*n+1
                  if(new_mutant_AMR_plasmids + plasmids_in_cell[2] > 0) col <- new_mutant_AMR_plasmids + plasmids_in_cell[2] else col <- 2*n+1
                  
                  if(!(0 %in% plasmids_in_cell)) {
                    probabilities_plasmids_after_duplication[row, col] <- 
                      ddirmnom(x = c(new_CRISPR_plasmids, new_mutant_AMR_plasmids, n - new_CRISPR_plasmids - new_mutant_AMR_plasmids - additional_CRISPR_plasmids), size = n- additional_CRISPR_plasmids, alpha = c(plasmids_in_cell[1],plasmids_in_cell[2] * (1 - competition_CRISPR_plasmid)* (1 - competition_AMR_plasmid * (CRISPR_type == "AMR ORI (mutations against AMR only)")), plasmids_in_cell[3] * (1 - competition_CRISPR_plasmid) * (1 - competition_AMR_plasmid)), log = FALSE)
                  } else if (plasmids_in_cell[1] == 0) {
                    if (plasmids_in_cell[2] > 0 && plasmids_in_cell[3] > 0)
                      probabilities_plasmids_after_duplication[2*n+1, col] <- 
                        ddirmnom(x = c(new_mutant_AMR_plasmids, n - new_mutant_AMR_plasmids - additional_CRISPR_plasmids), size = n- additional_CRISPR_plasmids, alpha = c(plasmids_in_cell[2] * (1 - competition_CRISPR_plasmid), plasmids_in_cell[3] * (1 - competition_CRISPR_plasmid)), log = FALSE)
                    else if (plasmids_in_cell[2] == 0)
                      probabilities_plasmids_after_duplication[2*n+1, 2*n+1] <- 1
                    else if (plasmids_in_cell[3] == 0)
                      probabilities_plasmids_after_duplication[2*n+1, 2*n] <- 1
                  } else if (plasmids_in_cell[1] > 0 && plasmids_in_cell[2] == 0) {
                    if (plasmids_in_cell[3] > 0)
                      probabilities_plasmids_after_duplication[row, 2*n+1] <- 
                        ddirmnom(x = c(new_CRISPR_plasmids, n - new_CRISPR_plasmids - additional_CRISPR_plasmids), size = n- additional_CRISPR_plasmids, alpha = c(plasmids_in_cell[1], plasmids_in_cell[3] * (1 - competition_CRISPR_plasmid) * (1 - competition_AMR_plasmid)), log = FALSE)
                    else if (plasmids_in_cell[3] == 0)
                      probabilities_plasmids_after_duplication[2*n, 2*n+1] <- 1
                  } else if (plasmids_in_cell[1] > 0 && plasmids_in_cell[2] > 0 && plasmids_in_cell[3] == 0)
                    probabilities_plasmids_after_duplication[row, col] <- 
                      ddirmnom(x = c(new_CRISPR_plasmids, new_mutant_AMR_plasmids), size = n- additional_CRISPR_plasmids, alpha = c(plasmids_in_cell[1],plasmids_in_cell[2] * (1 - competition_CRISPR_plasmid)* (1 - competition_AMR_plasmid * (CRISPR_type == "AMR ORI (mutations against AMR only)"))), log = FALSE)
                  
                }
              }
              
            } else if (type_plasmid_duplication == "regular") {
              
              # Each plasmid is replicated equally often. If this is not possible due to being more plasmids in the cell, each is replicated at most one time and chosen randomly (hypergeometrically).
              
              for (new_CRISPR_plasmids in 0:(n - additional_CRISPR_plasmids)) {
                for (new_mutant_AMR_plasmids in 0:(n-additional_CRISPR_plasmids)) {
                  if(new_CRISPR_plasmids + plasmids_in_cell[1] > 0) row <- new_CRISPR_plasmids + plasmids_in_cell[1] else row <- 2*n+1
                  if(new_mutant_AMR_plasmids + plasmids_in_cell[2] > 0) col <- new_mutant_AMR_plasmids + plasmids_in_cell[2] else col <- 2*n+1
                  
                  probabilities_plasmids_after_duplication[row, col] <- 
                    dmvhyper(x = c(new_CRISPR_plasmids, new_mutant_AMR_plasmids, n - new_CRISPR_plasmids - new_mutant_AMR_plasmids - additional_CRISPR_plasmids), n = plasmids_in_cell, k = n - additional_CRISPR_plasmids, log = FALSE)
                }
              }
            } else {print("Error (wrong type)")}
            
            # Calculation of the distribution of the number of plasmids to the daughter cells after plasmid duplication:
            
            probabilities_plasmids_in_daugther_cells$probability <- probabilities_plasmids_in_daugther_cells$probability + 
              probabilities_plasmids_after_duplication[cbind(probabilities_plasmids_in_daugther_cells$CRISPR_parent_after_doubling,probabilities_plasmids_in_daugther_cells$mutant_AMR_parent_after_doubling)] *
              dmvhyper(x = cbind(probabilities_plasmids_in_daugther_cells$CRISPR_first_daugther%%(n+1), probabilities_plasmids_in_daugther_cells$mutant_AMR_first_daugther%%(n+1), 
                                 n - probabilities_plasmids_in_daugther_cells$CRISPR_first_daugther%%(n+1) - probabilities_plasmids_in_daugther_cells$mutant_AMR_first_daugther%%(n+1)),
                       n = cbind(probabilities_plasmids_in_daugther_cells$CRISPR_parent_after_doubling%%(2*n+1), probabilities_plasmids_in_daugther_cells$mutant_AMR_parent_after_doubling%%(2*n+1), 
                                 2*n - probabilities_plasmids_in_daugther_cells$CRISPR_parent_after_doubling%%(2*n+1) - probabilities_plasmids_in_daugther_cells$mutant_AMR_parent_after_doubling%%(2*n+1)), k = n, log = FALSE) *
              ((additional_CRISPR_plasmids < n)*probability_u^(additional_CRISPR_plasmids)*(1-probability_u) + (additional_CRISPR_plasmids == n) * probability_u^n/(1-probability_u)*(1-probability_u))
            
            
          }
          
          # Ordering of the probabilities and removing non-occuring cases:
          
          probabilities_plasmids_in_daugther_cells <- probabilities_plasmids_in_daugther_cells[order(probabilities_plasmids_in_daugther_cells$probability,decreasing=TRUE),]
          
          probabilities_plasmids_in_daugther_cells <- subset(probabilities_plasmids_in_daugther_cells, probabilities_plasmids_in_daugther_cells$probability > 0)
          
          result[[i]][[mother_cell_CRISPR_plasmids]][[mother_cell_mutated_AMR_plasmids]] <- probabilities_plasmids_in_daugther_cells[,c(1:4,7)]
        }
        
      }
    }
    
    return(result)
    
  }
  
  else if (CRISPR_type == "Compatible Silencing" || CRISPR_type == "Compatible No Effect") {
    
    # Create a list for the result: The first entry in the list is the copy number, the second is the number of mutated AMR plasmids of the mother cell.
    
    result <- vector("list", length = length(observed_copy_numbers))
    
    for (l in 1:length(observed_copy_numbers)) {
      n <- observed_copy_numbers[l]
      result[[l]] <- vector("list", length = n + 1)
      
      for (i in 1:(n + 1)) {
        result[[l]][[i]] <- vector("list", length = n + 1)
      }
    }
    
    # Calculate the probabilities for the daughter cells for all mother cells and observed copy numbers:
    
    for (i in 1:length(observed_copy_numbers)) {
      
      n <- observed_copy_numbers[i]
      birth_death_parameter <- birth_death_parameters(n, lambda_max, CRISPR_type, dominance_function, lambda_min)
      
      for (mother_cell_mutated_AMR_plasmids in 1:(n+1)) {
        
        probability_u <- u / (birth_death_parameter[[1]][mother_cell_mutated_AMR_plasmids] + u)
        
        # Calculate the distribution of the number of plasmids after duplication.
        
        probabilities_plasmids_after_duplication <- rep(0, 2*n + 1)
        
        plasmids_in_cell <- c(mother_cell_mutated_AMR_plasmids%%(n+1), n - mother_cell_mutated_AMR_plasmids%%(n+1))
        
        if (type_plasmid_duplication == "random") {
          for (new_mutant_AMR_plasmids in 0:n) {
            if(new_mutant_AMR_plasmids + plasmids_in_cell[1] > 0) col <- new_mutant_AMR_plasmids + plasmids_in_cell[1] else col <- 2*n+1
            
            if(!(0 %in% plasmids_in_cell)) {
              probabilities_plasmids_after_duplication[col] <- 
                ddirmnom(x = c(new_mutant_AMR_plasmids, n - new_mutant_AMR_plasmids), size = n, alpha = c(plasmids_in_cell[1], plasmids_in_cell[2] * (1- competition_AMR_plasmid)), log = FALSE)
              
            } else if (plasmids_in_cell[1] == 0) { # Homozygote wildtype cells
              
              probabilities_plasmids_after_duplication[2*n+1] <- 1
              
            } else if (plasmids_in_cell[2] == 0) { # Homozygote mutant cells
              
              probabilities_plasmids_after_duplication[2*n] <- 1 
                
            }
            
          }
      } else if (type_plasmid_duplication == "regular") {
          
          # Each plasmid is duplicated exactly once.
          
          if (plasmids_in_cell[1] != 0) {
            probabilities_plasmids_after_duplication[2*plasmids_in_cell[1]] <- 1
          } else probabilities_plasmids_after_duplication[2*n+1] <- 1
            
          
      } else {print("Error (wrong type)")}
        
        # Now after the plasmids are duplicated, they are split up and given to the daughter cells:
        
        # First we create a matrix of all possible combinations of daughter cells:
        probabilities_plasmids_in_daugther_cells <- expand.grid(
          mutant_AMR_first_daugther = 1:(n+1),
          mutant_AMR_second_daugther = 1:(n+1)
        )
        
        probabilities_plasmids_in_daugther_cells$mutant_AMR_parent_after_doubling <- probabilities_plasmids_in_daugther_cells$mutant_AMR_first_daugther%%(n+1) + probabilities_plasmids_in_daugther_cells$mutant_AMR_second_daugther%%(n+1)
        
        probabilities_plasmids_in_daugther_cells$mutant_AMR_parent_after_doubling <- probabilities_plasmids_in_daugther_cells$mutant_AMR_parent_after_doubling + (2*n+1)*(probabilities_plasmids_in_daugther_cells$mutant_AMR_parent_after_doubling == 0)
        
        probabilities_plasmids_in_daugther_cells$probability_for_parent <- probabilities_plasmids_after_duplication[probabilities_plasmids_in_daugther_cells$mutant_AMR_parent_after_doubling]
        
        probabilities_plasmids_in_daugther_cells$probability <- probabilities_plasmids_in_daugther_cells$probability_for_parent *
          dmvhyper(x = cbind(probabilities_plasmids_in_daugther_cells$mutant_AMR_first_daugther%%(n+1), 
                             n - probabilities_plasmids_in_daugther_cells$mutant_AMR_first_daugther%%(n+1)),
                   n = cbind(probabilities_plasmids_in_daugther_cells$mutant_AMR_parent_after_doubling%%(2*n+1), 
                             2*n - probabilities_plasmids_in_daugther_cells$mutant_AMR_parent_after_doubling%%(2*n+1)), k = n, log = FALSE)
        
        # Sort by probability
        
        probabilities_plasmids_in_daugther_cells <- probabilities_plasmids_in_daugther_cells[order(probabilities_plasmids_in_daugther_cells$probability,decreasing=TRUE),]
        
        probabilities_plasmids_in_daugther_cells <- subset(probabilities_plasmids_in_daugther_cells, probabilities_plasmids_in_daugther_cells$probability > 0)
        
        # Save as result
        
        result[[i]][[mother_cell_mutated_AMR_plasmids]] <- probabilities_plasmids_in_daugther_cells[,c(1:2,5)]
      }
    }
    
    return(result)
    
  }
  
  else if (CRISPR_type == "Incompatible Cleaving") {
    
    # Create a list for the result: The first entry in the list is the copy number, the second is the number of mutated AMR plasmids of the mother cell.
    result <- vector("list", length = length(observed_copy_numbers))
    
    for (l in 1:length(observed_copy_numbers)) {
      n <- observed_copy_numbers[l]
      result[[l]] <- vector("list", length = n + 1)
      
      for (i in 1:(n + 1)) {
        result[[l]][[i]] <- vector("list", length = n + 1)
        
      }
    }
    
    # Calculate the probabilities for the daughter cells for all mother cells and observed copy numbers:
    
    for (i in 1:length(observed_copy_numbers)) {
      
      n <- observed_copy_numbers[i]
      birth_death_parameter <- birth_death_parameters(n, lambda_max, CRISPR_type, dominance_function, lambda_min)
    
      for (mother_cell_mutated_AMR_plasmids in 1:(n+1)) {
        
        probability_u <- u / (birth_death_parameter[[1]][mother_cell_mutated_AMR_plasmids] + u)
        
        # Calculate the distribution of the number of plasmids after duplication.
        
        probabilities_plasmids_after_duplication <- rep(0, 2*n + 1)
        
        plasmids_in_cell <- c(mother_cell_mutated_AMR_plasmids%%(n+1), n - mother_cell_mutated_AMR_plasmids%%(n+1))
        
        if (type_plasmid_duplication == "random") {
          for (new_mutant_AMR_plasmids in 0:n) {
            if(new_mutant_AMR_plasmids + plasmids_in_cell[1] > 0) col <- new_mutant_AMR_plasmids + plasmids_in_cell[1] else col <- 2*n+1
            
            if(!(0 %in% plasmids_in_cell)) {
              for (additional_CRISPR_plasmids in 0:n) {
                probabilities_plasmids_after_duplication[col] <- probabilities_plasmids_after_duplication[col] +
                  ddirmnom(x = c(new_mutant_AMR_plasmids, n - new_mutant_AMR_plasmids - additional_CRISPR_plasmids), size = n - additional_CRISPR_plasmids, alpha = c(plasmids_in_cell[1]* (1- competition_CRISPR_plasmid), (plasmids_in_cell[2] + additional_CRISPR_plasmids)), log = FALSE) * 
                  ((additional_CRISPR_plasmids < n)*probability_u^(additional_CRISPR_plasmids)*(1-probability_u) + (additional_CRISPR_plasmids == n) * probability_u^n/(1-probability_u)*(1-probability_u))
              }
              
            } else if (plasmids_in_cell[1] == 0) { # Homozygote wildtype cells
              
              probabilities_plasmids_after_duplication[2*n+1] <- 1
              
            } else if (plasmids_in_cell[2] == 0) { # Homozygote mutant cells

              for (additional_CRISPR_plasmids in 0:n) {
                if (additional_CRISPR_plasmids == 0 && new_mutant_AMR_plasmids == n) {
                  probabilities_plasmids_after_duplication[col] <- probabilities_plasmids_after_duplication[col] +
                    ((additional_CRISPR_plasmids < n)*probability_u^(additional_CRISPR_plasmids)*(1-probability_u) + (additional_CRISPR_plasmids == n) * probability_u^n/(1-probability_u)*(1-probability_u))  
                } else if (additional_CRISPR_plasmids > 0){
                  probabilities_plasmids_after_duplication[col] <- probabilities_plasmids_after_duplication[col] +
                    ddirmnom(x = c(new_mutant_AMR_plasmids, n - new_mutant_AMR_plasmids - additional_CRISPR_plasmids), size = n - additional_CRISPR_plasmids, alpha = c(plasmids_in_cell[1]* (1- competition_CRISPR_plasmid), (plasmids_in_cell[2] + additional_CRISPR_plasmids)), log = FALSE) * 
                    ((additional_CRISPR_plasmids < n)*probability_u^(additional_CRISPR_plasmids)*(1-probability_u) + (additional_CRISPR_plasmids == n) * probability_u^n/(1-probability_u)*(1-probability_u))
                }
              
              }
            }
          }
        } else if (type_plasmid_duplication == "regular") {
          
          # Each plasmid is duplicated exactly once.
          
         
          for (additional_CRISPR_plasmids in 0:n) {
            for (new_mutant_AMR_plasmids in 0:(n-additional_CRISPR_plasmids)) {
              if(new_mutant_AMR_plasmids + plasmids_in_cell[1] > 0) col <- new_mutant_AMR_plasmids + plasmids_in_cell[1] else col <- 2*n+1
              
              probabilities_plasmids_after_duplication[col] <- probabilities_plasmids_after_duplication[col] + 
                dmvhyper(x = c(new_mutant_AMR_plasmids, n - new_mutant_AMR_plasmids - additional_CRISPR_plasmids), n = plasmids_in_cell, k = n - additional_CRISPR_plasmids, log = FALSE) *
                ((additional_CRISPR_plasmids < n)*probability_u^(additional_CRISPR_plasmids)*(1-probability_u) + (additional_CRISPR_plasmids == n) * probability_u^n/(1-probability_u)*(1-probability_u))
            }
          
          if (plasmids_in_cell[1] != 0) probabilities_plasmids_after_duplication[2*plasmids_in_cell[1]] <- 1
          else probabilities_plasmids_after_duplication[2*n+1] <- 1
            
          }
          
        } else {print("Error (wrong type)")}
        
        # Now after the plasmids are duplicated, they are split up and given to the daughter cells:
        
        # First we create a matrix of all possible combinations of daughter cells:
        probabilities_plasmids_in_daugther_cells <- expand.grid(
          mutant_AMR_first_daugther = 1:(n+1),
          mutant_AMR_second_daugther = 1:(n+1)
        )
        
        probabilities_plasmids_in_daugther_cells$mutant_AMR_parent_after_doubling <- probabilities_plasmids_in_daugther_cells$mutant_AMR_first_daugther%%(n+1) + probabilities_plasmids_in_daugther_cells$mutant_AMR_second_daugther%%(n+1)
        
        probabilities_plasmids_in_daugther_cells$mutant_AMR_parent_after_doubling <- probabilities_plasmids_in_daugther_cells$mutant_AMR_parent_after_doubling + (2*n+1)*(probabilities_plasmids_in_daugther_cells$mutant_AMR_parent_after_doubling == 0)
        
        probabilities_plasmids_in_daugther_cells$probability_for_parent <- probabilities_plasmids_after_duplication[probabilities_plasmids_in_daugther_cells$mutant_AMR_parent_after_doubling]
        
        probabilities_plasmids_in_daugther_cells$probability <- probabilities_plasmids_in_daugther_cells$probability_for_parent *
          dmvhyper(x = cbind(probabilities_plasmids_in_daugther_cells$mutant_AMR_first_daugther%%(n+1), 
                             n - probabilities_plasmids_in_daugther_cells$mutant_AMR_first_daugther%%(n+1)),
                   n = cbind(probabilities_plasmids_in_daugther_cells$mutant_AMR_parent_after_doubling%%(2*n+1), 
                             2*n - probabilities_plasmids_in_daugther_cells$mutant_AMR_parent_after_doubling%%(2*n+1)), k = n, log = FALSE)
        
        # Sort by probability
        
        probabilities_plasmids_in_daugther_cells <- probabilities_plasmids_in_daugther_cells[order(probabilities_plasmids_in_daugther_cells$probability,decreasing=TRUE),]
        
        probabilities_plasmids_in_daugther_cells <- subset(probabilities_plasmids_in_daugther_cells, probabilities_plasmids_in_daugther_cells$probability > 0)
        
        # Save as result
        
        result[[i]][[mother_cell_mutated_AMR_plasmids]] <- probabilities_plasmids_in_daugther_cells[,c(1:2,5)]
        
      }
    }
    
    
    return(result)
    
    
  }
  
  else if (CRISPR_type == "Compatible Cleaving") {
    
    
    # Create a list for the result: The first entry in the list is the copy number, the second is the number of mutated AMR plasmids of the mother cell.
    result <- vector("list", length = length(observed_copy_numbers))
    
    for (l in 1:length(observed_copy_numbers)) {
      n <- observed_copy_numbers[l]
      result[[l]] <- vector("list", length = n + 1)
      
      for (i in 1:(n + 1)) {
        result[[l]][[i]] <- vector("list", length = n + 1)
        
      }
    }
    # The first entry in the list is the copy number, the second is the number of mutated AMR plasmids of the mother cell. 
    
    for (i in 1:length(observed_copy_numbers)) {
      
      n <- observed_copy_numbers[i]
      
      for (mother_cell_mutated_AMR_plasmids in 1:(n+1)) {
        
        probability_u <- 0
        
        # Calculate the distribution of the number of plasmids after duplication.
        
        probabilities_plasmids_after_duplication <- rep(0, 2*n + 1)
        
        plasmids_in_cell <- c(mother_cell_mutated_AMR_plasmids%%(n+1), n - mother_cell_mutated_AMR_plasmids%%(n+1))
        
        if (type_plasmid_duplication == "random") {
          for (new_mutant_AMR_plasmids in 0:n) {
            if(new_mutant_AMR_plasmids + plasmids_in_cell[1] > 0) col <- new_mutant_AMR_plasmids + plasmids_in_cell[1] else col <- 2*n+1
            
            if(!(0 %in% plasmids_in_cell)) {
              for (additional_CRISPR_plasmids in 0:n) {
                probabilities_plasmids_after_duplication[col] <- probabilities_plasmids_after_duplication[col] +
                  ddirmnom(x = c(new_mutant_AMR_plasmids, n - new_mutant_AMR_plasmids - additional_CRISPR_plasmids), size = n - additional_CRISPR_plasmids, alpha = c(plasmids_in_cell[1]* (1- competition_CRISPR_plasmid), (plasmids_in_cell[2] + additional_CRISPR_plasmids)), log = FALSE) * 
                  ((additional_CRISPR_plasmids < n)*probability_u^(additional_CRISPR_plasmids)*(1-probability_u) + (additional_CRISPR_plasmids == n) * probability_u^n/(1-probability_u)*(1-probability_u))
              }
              
            } else if (plasmids_in_cell[1] == 0) { # Homozygote wildtype cells
              
              probabilities_plasmids_after_duplication[2*n+1] <- 1
              
            } else if (plasmids_in_cell[2] == 0) { # Homozygote mutant cells
              
              for (additional_CRISPR_plasmids in 0:n) {
                if (additional_CRISPR_plasmids == 0 && new_mutant_AMR_plasmids == n) {
                  probabilities_plasmids_after_duplication[col] <- probabilities_plasmids_after_duplication[col] +
                    ((additional_CRISPR_plasmids < n)*probability_u^(additional_CRISPR_plasmids)*(1-probability_u) + (additional_CRISPR_plasmids == n) * probability_u^n/(1-probability_u)*(1-probability_u))  
                } else if (additional_CRISPR_plasmids > 0){
                  probabilities_plasmids_after_duplication[col] <- probabilities_plasmids_after_duplication[col] +
                    ddirmnom(x = c(new_mutant_AMR_plasmids, n - new_mutant_AMR_plasmids - additional_CRISPR_plasmids), size = n - additional_CRISPR_plasmids, alpha = c(plasmids_in_cell[1]* (1- competition_CRISPR_plasmid), (plasmids_in_cell[2] + additional_CRISPR_plasmids)), log = FALSE) * 
                    ((additional_CRISPR_plasmids < n)*probability_u^(additional_CRISPR_plasmids)*(1-probability_u) + (additional_CRISPR_plasmids == n) * probability_u^n/(1-probability_u)*(1-probability_u))
                }
                
              }
            }
          }
        } else if (type_plasmid_duplication == "regular") {
          
          # Each plasmid is duplicated exactly once.
          
          
          for (additional_CRISPR_plasmids in 0:n) {
            for (new_mutant_AMR_plasmids in 0:(n-additional_CRISPR_plasmids)) {
              if(new_mutant_AMR_plasmids + plasmids_in_cell[1] > 0) col <- new_mutant_AMR_plasmids + plasmids_in_cell[1] else col <- 2*n+1
              
              probabilities_plasmids_after_duplication[col] <- probabilities_plasmids_after_duplication[col] + 
                dmvhyper(x = c(new_mutant_AMR_plasmids, n - new_mutant_AMR_plasmids - additional_CRISPR_plasmids), n = plasmids_in_cell, k = n - additional_CRISPR_plasmids, log = FALSE) *
                ((additional_CRISPR_plasmids < n)*probability_u^(additional_CRISPR_plasmids)*(1-probability_u) + (additional_CRISPR_plasmids == n) * probability_u^n/(1-probability_u)*(1-probability_u))
            }
            
            
            
            if (plasmids_in_cell[1] != 0) probabilities_plasmids_after_duplication[2*plasmids_in_cell[1]] <- 1
            else probabilities_plasmids_after_duplication[2*n+1] <- 1
            
          }
          
        } else {print("Error (wrong type)")}
        
        # Now after the plasmids are duplicated, they are split up and given to the daughter cells:
        
        # First we create a matrix of all possible combinations of daughter cells:
        probabilities_plasmids_in_daugther_cells <- expand.grid(
          mutant_AMR_first_daugther = 1:(n+1),
          mutant_AMR_second_daugther = 1:(n+1)
        )
        
        probabilities_plasmids_in_daugther_cells$mutant_AMR_parent_after_doubling <- probabilities_plasmids_in_daugther_cells$mutant_AMR_first_daugther%%(n+1) + probabilities_plasmids_in_daugther_cells$mutant_AMR_second_daugther%%(n+1)
        
        probabilities_plasmids_in_daugther_cells$mutant_AMR_parent_after_doubling <- probabilities_plasmids_in_daugther_cells$mutant_AMR_parent_after_doubling + (2*n+1)*(probabilities_plasmids_in_daugther_cells$mutant_AMR_parent_after_doubling == 0)
        
        probabilities_plasmids_in_daugther_cells$probability_for_parent <- probabilities_plasmids_after_duplication[probabilities_plasmids_in_daugther_cells$mutant_AMR_parent_after_doubling]
        
        probabilities_plasmids_in_daugther_cells$probability <- probabilities_plasmids_in_daugther_cells$probability_for_parent *
          dmvhyper(x = cbind(probabilities_plasmids_in_daugther_cells$mutant_AMR_first_daugther%%(n+1), 
                             n - probabilities_plasmids_in_daugther_cells$mutant_AMR_first_daugther%%(n+1)),
                   n = cbind(probabilities_plasmids_in_daugther_cells$mutant_AMR_parent_after_doubling%%(2*n+1), 
                             2*n - probabilities_plasmids_in_daugther_cells$mutant_AMR_parent_after_doubling%%(2*n+1)), k = n, log = FALSE)
        
        # Sort by probability
        
        probabilities_plasmids_in_daugther_cells <- probabilities_plasmids_in_daugther_cells[order(probabilities_plasmids_in_daugther_cells$probability,decreasing=TRUE),]
        
        probabilities_plasmids_in_daugther_cells <- subset(probabilities_plasmids_in_daugther_cells, probabilities_plasmids_in_daugther_cells$probability > 0)
        
        # Save as result
        
        result[[i]][[mother_cell_mutated_AMR_plasmids]] <- probabilities_plasmids_in_daugther_cells[,c(1:2,5)]
        
      }
    }
    
    return(result)

  }
  
  else {print("Error (wrong type)")}
  
}

