
# Calculation of the probabilities of a cell with number_of_mutations mutated AMR plasmids, after it gained a CRISPR encoding plasmid.

# To avoid that there are more than n (copy number) plasmids in the cell, we wait for one generation to pass and calculate the probabilities of all cell types, which the two daughter cells can have. 
# Otherwise an incompatible AMR plasmid could lead to n+1 plasmids in total. The results for Fig 2 in the paper will be calculated differently.

# OUTPUT: Probability of each cell type after CRISPR introduction.


first_generation_probabilities <- function(number_of_mutations, n, CRISPR_type, type_plasmid_duplication, s_max, s_0, u, dominance_function, competition_CRISPR_plasmid, competition_ABR_plasmid) {
  
  if (CRISPR_type == "Incompatible Silencing" || CRISPR_type == "Incompatible No Effect") {
    
    # Number of plasmids in the cell:
    
    plasmids_in_cell <- c(1, number_of_mutations%%(n+1), n - 1 - number_of_mutations%%(n+1) + 1)
    
    # Birth parameters for the moment of CRISPR introduction:
    
    lambda <- rep(1 + s_max, n + 1)
    mu <- rep(1, n+1)
    
    if (CRISPR_type == "Incompatible Silencing") {
      
      if (dominance_function == "dominant") {
        
        lambda[n+1] <- 1 - s_0
        
      } else if (dominance_function == "recessive") {
        
        lambda <- rep(1 - s_0, n + 1)
        
        lambda[n] <- 1 + s_max
        
      } else if (dominance_function == "linear") {
        
        lambda <- seq(from = 1 - s_0, to = 1 + s_max , length.out = n + 1) 
        
        lambda <- c(lambda[2:(n+1)], lambda[1])
        
      }
    }
    
    # Probability to get additional CRISPR encoding plasmids after the first one:  
    
    if (number_of_mutations > 0) {
      probability_u <- u / (u + lambda[number_of_mutations])
    } else {
      probability_u <- u / (u + lambda[n+1])
    }
    
    
    # Probabilities of the amount of plasmids after plasmid duplication (before cell division):
    
    probabilities_plasmids_in_daugther_cells <- expand.grid(
      CRISPR_first_daugther = 1:(n+1),
      mutant_ABR_first_daugther = 1:(n+1),
      CRISPR_second_daugther = 1:(n+1),
      mutant_ABR_second_daugther = 1:(n+1)
    )
    
    probabilities_plasmids_in_daugther_cells <- subset(probabilities_plasmids_in_daugther_cells, probabilities_plasmids_in_daugther_cells$CRISPR_first_daugther%%(n+1) + probabilities_plasmids_in_daugther_cells$mutant_ABR_first_daugther%%(n+1)  <= n &
                                                         probabilities_plasmids_in_daugther_cells$CRISPR_second_daugther%%(n+1) + probabilities_plasmids_in_daugther_cells$mutant_ABR_second_daugther%%(n+1)  <= n)
    
    probabilities_plasmids_in_daugther_cells$CRISPR_parent_after_doubling <- probabilities_plasmids_in_daugther_cells$CRISPR_first_daugther%%(n+1) + probabilities_plasmids_in_daugther_cells$CRISPR_second_daugther%%(n+1)
    probabilities_plasmids_in_daugther_cells$mutant_ABR_parent_after_doubling <- probabilities_plasmids_in_daugther_cells$mutant_ABR_first_daugther%%(n+1) + probabilities_plasmids_in_daugther_cells$mutant_ABR_second_daugther%%(n+1)
    
    probabilities_plasmids_in_daugther_cells$CRISPR_parent_after_doubling <- probabilities_plasmids_in_daugther_cells$CRISPR_parent_after_doubling + (2*n+1)*(probabilities_plasmids_in_daugther_cells$CRISPR_parent_after_doubling == 0)
    probabilities_plasmids_in_daugther_cells$mutant_ABR_parent_after_doubling <- probabilities_plasmids_in_daugther_cells$mutant_ABR_parent_after_doubling + (2*n+1)*(probabilities_plasmids_in_daugther_cells$mutant_ABR_parent_after_doubling == 0)
    
    probabilities_plasmids_in_daugther_cells$probability <- 0
    
    
    for (additional_CRISPR_plasmids in 1:n) {
      
      probabilities_plasmids_after_duplication <- matrix(data = 0, nrow = 2*n+1, ncol = 2*n+1)
      
      plasmids_in_cell <- c(additional_CRISPR_plasmids, number_of_mutations%%(n+1), 
                            n - number_of_mutations%%(n+1))
      
      
      if (type_plasmid_duplication == "random") {
        
        # Polya's urn model: Only one plasmid is replicated at a time. The plasmid is chosen randomly. The duplicated plasmid is added to the pool.
        
        for (new_CRISPR_plasmids in 0:(n-additional_CRISPR_plasmids)) {
          for (new_mutant_ABR_plasmids in 0:(n-additional_CRISPR_plasmids)) {
            if(new_CRISPR_plasmids + plasmids_in_cell[1] > 0) row <- new_CRISPR_plasmids + plasmids_in_cell[1] else row <- 2*n+1
            if(new_mutant_ABR_plasmids + plasmids_in_cell[2] > 0) col <- new_mutant_ABR_plasmids + plasmids_in_cell[2] else col <- 2*n+1
            
            if(!(0 %in% plasmids_in_cell)) {
              probabilities_plasmids_after_duplication[row, col] <- 
                ddirmnom(x = c(new_CRISPR_plasmids, new_mutant_ABR_plasmids, n - new_CRISPR_plasmids - new_mutant_ABR_plasmids - additional_CRISPR_plasmids), size = n- additional_CRISPR_plasmids, alpha = c(plasmids_in_cell[1],plasmids_in_cell[2] * (1 - competition_CRISPR_plasmid)* (1 - competition_ABR_plasmid * (CRISPR_type == "ABR ORI (mutations against ABR only)")), plasmids_in_cell[3] * (1 - competition_CRISPR_plasmid) * (1 - competition_ABR_plasmid)), log = FALSE)
            } else if (plasmids_in_cell[1] == 0) {
              if (plasmids_in_cell[2] > 0 && plasmids_in_cell[3] > 0)
                probabilities_plasmids_after_duplication[2*n+1, col] <- 
                  ddirmnom(x = c(new_mutant_ABR_plasmids, n - new_mutant_ABR_plasmids - additional_CRISPR_plasmids), size = n- additional_CRISPR_plasmids, alpha = c(plasmids_in_cell[2] * (1 - competition_CRISPR_plasmid), plasmids_in_cell[3] * (1 - competition_CRISPR_plasmid)), log = FALSE)
              else if (plasmids_in_cell[2] == 0)
                probabilities_plasmids_after_duplication[2*n+1, 2*n+1] <- 1
              else if (plasmids_in_cell[3] == 0)
                probabilities_plasmids_after_duplication[2*n+1, 2*n] <- 1
            } else if (plasmids_in_cell[1] > 0 && plasmids_in_cell[2] == 0) {
              if (plasmids_in_cell[3] > 0)
                probabilities_plasmids_after_duplication[row, 2*n+1] <- 
                  ddirmnom(x = c(new_CRISPR_plasmids, n - new_CRISPR_plasmids - additional_CRISPR_plasmids), size = n- additional_CRISPR_plasmids, alpha = c(plasmids_in_cell[1], plasmids_in_cell[3] * (1 - competition_CRISPR_plasmid) * (1 - competition_ABR_plasmid)), log = FALSE)
              else if (plasmids_in_cell[3] == 0)
                probabilities_plasmids_after_duplication[2*n, 2*n+1] <- 1
            } else if (plasmids_in_cell[1] > 0 && plasmids_in_cell[2] > 0 && plasmids_in_cell[3] == 0)
              probabilities_plasmids_after_duplication[row, col] <- 
                ddirmnom(x = c(new_CRISPR_plasmids, new_mutant_ABR_plasmids), size = n- additional_CRISPR_plasmids, alpha = c(plasmids_in_cell[1],plasmids_in_cell[2] * (1 - competition_CRISPR_plasmid))* (1 - competition_ABR_plasmid * (CRISPR_type == "ABR ORI (mutations against ABR only)")), log = FALSE)
            
          }
        }
        
      } else if (type_plasmid_duplication == "regular") {
        
        # Each plasmid is replicated equally often. If this is not possible due to being more plasmids in the cell, each is replicated at most one time and chosen randomly (hypergeometrically).
        
        for (new_CRISPR_plasmids in 0:(n - additional_CRISPR_plasmids)) {
          for (new_mutant_ABR_plasmids in 0:(n-additional_CRISPR_plasmids)) {
            if(new_CRISPR_plasmids + plasmids_in_cell[1] > 0) row <- new_CRISPR_plasmids + plasmids_in_cell[1] else row <- 2*n+1
            if(new_mutant_ABR_plasmids + plasmids_in_cell[2] > 0) col <- new_mutant_ABR_plasmids + plasmids_in_cell[2] else col <- 2*n+1
            
            probabilities_plasmids_after_duplication[row, col] <- 
              dmvhyper(x = c(new_CRISPR_plasmids, new_mutant_ABR_plasmids, n - new_CRISPR_plasmids - new_mutant_ABR_plasmids - additional_CRISPR_plasmids), n = plasmids_in_cell, k = n - additional_CRISPR_plasmids, log = FALSE)
          }
        }
      } else {print("Error (wrong type)")}
      
      
      # Calculation of the distribution of the number of plasmids to the daughter cells after plasmid duplication:
      
      probabilities_plasmids_in_daugther_cells$probability <- probabilities_plasmids_in_daugther_cells$probability + 
        probabilities_plasmids_after_duplication[cbind(probabilities_plasmids_in_daugther_cells$CRISPR_parent_after_doubling,probabilities_plasmids_in_daugther_cells$mutant_ABR_parent_after_doubling)] *
        dmvhyper(x = cbind(probabilities_plasmids_in_daugther_cells$CRISPR_first_daugther%%(n+1), probabilities_plasmids_in_daugther_cells$mutant_ABR_first_daugther%%(n+1), 
                           n - probabilities_plasmids_in_daugther_cells$CRISPR_first_daugther%%(n+1) - probabilities_plasmids_in_daugther_cells$mutant_ABR_first_daugther%%(n+1)),
                 n = cbind(probabilities_plasmids_in_daugther_cells$CRISPR_parent_after_doubling%%(2*n+1), probabilities_plasmids_in_daugther_cells$mutant_ABR_parent_after_doubling%%(2*n+1), 
                           2*n - probabilities_plasmids_in_daugther_cells$CRISPR_parent_after_doubling%%(2*n+1) - probabilities_plasmids_in_daugther_cells$mutant_ABR_parent_after_doubling%%(2*n+1)), k = n, log = FALSE) *
        ((additional_CRISPR_plasmids < n)*probability_u^(additional_CRISPR_plasmids-1)*(1-probability_u) + (additional_CRISPR_plasmids == n) * probability_u^(n-1)/(1-probability_u)*(1-probability_u))
      
      
    }
    
    
    probabilities_plasmids_in_daugther_cells <- probabilities_plasmids_in_daugther_cells[,c(1:4,7)]
    
    
    # Addition of the probability, that the cell dies out instead of replicating:
    
    plasmids_in_cell <- c(1, number_of_mutations, n - 1 - number_of_mutations + 1)
    
    if (number_of_mutations > 0) {
      probability_birth <- lambda[number_of_mutations] / (mu[number_of_mutations] + lambda[number_of_mutations])
    } else {
      probability_birth <- min(lambda) / (max(mu) + min(lambda))
    }
    
    probabilities_plasmids_in_daugther_cells$probability <- probabilities_plasmids_in_daugther_cells$probability * probability_birth 
    
    probabilities_plasmids_in_daugther_cells <- rbind(probabilities_plasmids_in_daugther_cells, c(rep(0, 4), 1- probability_birth))
    
    
    # Ordering of probabilities and removal of non-occurring cases:
    
    probabilities_plasmids_in_daugther_cells <- probabilities_plasmids_in_daugther_cells[order(probabilities_plasmids_in_daugther_cells$probability,decreasing=TRUE),]
    
    probabilities_plasmids_in_daugther_cells <- subset(probabilities_plasmids_in_daugther_cells, probabilities_plasmids_in_daugther_cells$probability > 0)
    
  } 
  
  else if (CRISPR_type == "Compatible Silencing" || CRISPR_type == "Incompatible Cleaving" || CRISPR_type == "Compatible Cleaving" || CRISPR_type == "Compatible No Effect") {
    
    # Calculation of the birth and death parameters for the cell with number_of_mutations mutated AMR plasmids:
    
    birth_death_parameter <- birth_death_parameters(n, s_max, CRISPR_type, dominance_function, s_0)
    
    # Calculation of the number of each plasmid type in the cell:
    
    if (CRISPR_type == "Incompatible Cleaving") { 
      
      plasmids_in_cell <- c(number_of_mutations%%(n+1), 1)
      
      if (type_plasmid_duplication == "regular") {while( sum(plasmids_in_cell) <= n) plasmids_in_cell <- plasmids_in_cell*2}
      
    } else if (CRISPR_type == "Compatible Cleaving") {
      
      plasmids_in_cell <- (number_of_mutations%%(n+1) != 0) * c(n, 0) + (number_of_mutations%%(n+1) == 0) * c(0, n)
    
    } else if (CRISPR_type == "Compatible Silencing") {
      
      plasmids_in_cell <- c(number_of_mutations%%(n+1), n - number_of_mutations%%(n+1))
    
    } else if (CRISPR_type == "Compatible No Effect") {
      
      plasmids_in_cell <- c(number_of_mutations%%(n+1), n - number_of_mutations%%(n+1))
      
    }
    
    # Probability of getting more than one CRISPR encoding plasmid:
    
    if (number_of_mutations > 0) {
      probability_u <- u / (u + birth_death_parameter[[1]][number_of_mutations])
    } else {
      probability_u <- u / (u + birth_death_parameter[[1]][n+1])
    }
    
    # Calculation of the amount of plasmids after the plasmids are duplicated for the next cell division:
    
    probabilities_plasmids_after_duplication <- rep(0, 2*n + 1)
    
    if (type_plasmid_duplication == "random") {
      
      # Polya's urn model: Only one plasmid is replicated at a time. The plasmid is chosen randomly. The duplicated plasmid is added to the pool.
      
      for (new_mutant_ABR_plasmids in 0:(2*n - plasmids_in_cell[1] - plasmids_in_cell[2])) {
        if(new_mutant_ABR_plasmids + plasmids_in_cell[1] > 0) col <- new_mutant_ABR_plasmids + plasmids_in_cell[1] else col <- 2*n + 1
        
        if(!(0 %in% plasmids_in_cell)) {
          if (CRISPR_type == "Compatible Silencing" || CRISPR_type == "Compatible No Effect") probabilities_plasmids_after_duplication[col] <- 
              ddirmnom(x = c(new_mutant_ABR_plasmids, (2*n - plasmids_in_cell[1] - plasmids_in_cell[2]) - new_mutant_ABR_plasmids), size = 2*n - plasmids_in_cell[1] - plasmids_in_cell[2], alpha = c(plasmids_in_cell[1], plasmids_in_cell[2] * (1- competition_ABR_plasmid)), log = FALSE)
          else if (CRISPR_type == "Incompatible Cleaving") {
            for (additional_CRISPR_plasmids in 0:(2*n - plasmids_in_cell[1] - plasmids_in_cell[2] - new_mutant_ABR_plasmids)) {
              probabilities_plasmids_after_duplication[col] <- probabilities_plasmids_after_duplication[col] +
                ddirmnom(x = c(new_mutant_ABR_plasmids, 2*n - plasmids_in_cell[1] - plasmids_in_cell[2]- new_mutant_ABR_plasmids - additional_CRISPR_plasmids), size = 2*n - plasmids_in_cell[1] - plasmids_in_cell[2] - additional_CRISPR_plasmids, alpha = c(plasmids_in_cell[1]* (1- competition_CRISPR_plasmid), (plasmids_in_cell[2] + additional_CRISPR_plasmids) ), log = FALSE) * 
                ((additional_CRISPR_plasmids < (2*n - plasmids_in_cell[1] - plasmids_in_cell[2]))*probability_u^(additional_CRISPR_plasmids)*(1-probability_u) + (additional_CRISPR_plasmids == (2*n - plasmids_in_cell[1] - plasmids_in_cell[2])) * probability_u^(2*n - plasmids_in_cell[1] - plasmids_in_cell[2])/(1-probability_u)*(1-probability_u))
            } 
          }
          
        } else if (plasmids_in_cell[1] == 0) {
          probabilities_plasmids_after_duplication[2*n+1] <- 1 
        } else if (plasmids_in_cell[2] == 0) {
          probabilities_plasmids_after_duplication[2*n] <- 1
        }
      }
    } else if (type_plasmid_duplication == "regular") {
      
      # Each plasmid is replicated equally often. If this is not possible due to being more plasmids in the cell, each is replicated at most one time and chosen randomly (hypergeometrically).
      
      if (CRISPR_type == "Compatible Silencing" || CRISPR_type == "Compatible No Effect") {
        
        if (plasmids_in_cell[1] != 0) probabilities_plasmids_after_duplication[2*plasmids_in_cell[1]] <- 1
        else probabilities_plasmids_after_duplication[2*n+1] <- 1
        
      } else if (CRISPR_type == "Incompatible Cleaving" || CRISPR_type == "Compatible Cleaving" || CRISPR_type == "Compatible No Effect") {
        
        for (new_mutant_ABR_plasmids in 0:(2*n-sum(plasmids_in_cell))) {
          if(new_mutant_ABR_plasmids + plasmids_in_cell[1] > 0) col <- new_mutant_ABR_plasmids + plasmids_in_cell[1] else col <- 2*n+1
          
          # If not every plasmid can be copied once, the plasmids are chosen randomly but at most once.
          
          probabilities_plasmids_after_duplication[col] <- 
            dmvhyper(x = c(new_mutant_ABR_plasmids, 2*n - new_mutant_ABR_plasmids - sum(plasmids_in_cell)), n = plasmids_in_cell, k = 2*n-sum(plasmids_in_cell), log = FALSE)
        }
      } 
    
  } else {print("Error (wrong type)")}
  
  # Calculation of the distribution of the number of plasmids to the daughter cells after plasmid duplication:
    
  probabilities_plasmids_in_daugther_cells <- expand.grid(
    mutant_ABR_first_daugther = 1:(n+1),
    mutant_ABR_second_daugther = 1:(n+1)
  )
  
  probabilities_plasmids_in_daugther_cells$mutant_ABR_parent_after_doubling <- probabilities_plasmids_in_daugther_cells$mutant_ABR_first_daugther%%(n+1) + probabilities_plasmids_in_daugther_cells$mutant_ABR_second_daugther%%(n+1)
  
  probabilities_plasmids_in_daugther_cells$mutant_ABR_parent_after_doubling <- probabilities_plasmids_in_daugther_cells$mutant_ABR_parent_after_doubling + (2*n+1)*(probabilities_plasmids_in_daugther_cells$mutant_ABR_parent_after_doubling == 0)
  
  probabilities_plasmids_in_daugther_cells$probability_for_parent <- probabilities_plasmids_after_duplication[probabilities_plasmids_in_daugther_cells$mutant_ABR_parent_after_doubling]
  
  probabilities_plasmids_in_daugther_cells$probability <- probabilities_plasmids_in_daugther_cells$probability_for_parent *
    dmvhyper(x = cbind(probabilities_plasmids_in_daugther_cells$mutant_ABR_first_daugther%%(n+1), 
                       n - probabilities_plasmids_in_daugther_cells$mutant_ABR_first_daugther%%(n+1)),
             n = cbind(probabilities_plasmids_in_daugther_cells$mutant_ABR_parent_after_doubling%%(2*n+1), 
                       2*n - probabilities_plasmids_in_daugther_cells$mutant_ABR_parent_after_doubling%%(2*n+1)), k = n, log = FALSE)
  
  probabilities_plasmids_in_daugther_cells <- probabilities_plasmids_in_daugther_cells[,c(1:2,5)] # Reduction of the output
  
  
  # Addition of the probability, that the cell dies out instead of replicating:
  
  if (CRISPR_type == "Compatible Silencing") {
    
    if (number_of_mutations > 0) {
      probability_birth <- birth_death_parameter[[1]][number_of_mutations] / (birth_death_parameter[[2]][number_of_mutations] + birth_death_parameter[[1]][number_of_mutations])
    } else {
      probability_birth <- min(birth_death_parameter[[1]]) / (max(birth_death_parameter[[2]]) + min(birth_death_parameter[[1]]))
    }
    
  } else if (CRISPR_type == "Incompatible Cleaving") {
    proportion_ABR <- plasmids_in_cell[1]/(plasmids_in_cell[1] + plasmids_in_cell[2])
    
    probability_birth <- (1-proportion_ABR%%1) * birth_death_parameter[[1]][max(1,floor(proportion_ABR*n))] / (birth_death_parameter[[2]][max(1,floor(proportion_ABR*n))] + birth_death_parameter[[1]][max(1,floor(proportion_ABR*n))]) +
      (proportion_ABR%%1) * birth_death_parameter[[1]][max(1,ceiling(proportion_ABR*n))] / (birth_death_parameter[[2]][max(1,ceiling(proportion_ABR*n))] + birth_death_parameter[[1]][max(1,ceiling(proportion_ABR*n))]) # Average result
    
  } else if (CRISPR_type == "Compatible Cleaving" || CRISPR_type == "Compatible No Effect") {
    probability_birth <- birth_death_parameter[[1]][n] / (birth_death_parameter[[2]][n] + birth_death_parameter[[1]][n]) 
  }
  
  probabilities_plasmids_in_daugther_cells$probability <- probabilities_plasmids_in_daugther_cells$probability * probability_birth 
  
  probabilities_plasmids_in_daugther_cells <- rbind(probabilities_plasmids_in_daugther_cells, c(rep(0, 2), 1- probability_birth))
  
  
  # Ordering of the probabilities and removing non-occuring cases:
  
  probabilities_plasmids_in_daugther_cells <- probabilities_plasmids_in_daugther_cells[order(probabilities_plasmids_in_daugther_cells$probability,decreasing=TRUE),]
  
  probabilities_plasmids_in_daugther_cells <- subset(probabilities_plasmids_in_daugther_cells, probabilities_plasmids_in_daugther_cells$probability > 0)
  
}
  
  return(probabilities_plasmids_in_daugther_cells)

}
