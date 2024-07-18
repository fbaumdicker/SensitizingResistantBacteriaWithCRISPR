# **Sensitizing Resistant Bacteria With CRISPR**

## **Required Software:**
- R
- Python with packages numpy and scipy (can be installed using the file `Standing_genetic_variation.Rmd`)

Before the R code can be run, it is necessary to run every section of `Standing_genetic_variation.Rmd` beforehand. This needs to be repeated if R is restarted.

## **Summary of All Files:**

- **`Birth-death-parameters.R`**  
  Calculation of the birth and death parameters in Phase 3.

- **`Cell_division_plasmid_probabilities.R`**  
  Calculation of the probability of all daughter cell types, given the plasmid composition in the mother cell. De novo mutations (after CRISPR introduction) are not included.

- **`Figure2.R`**  
  Calculation of the data needed for Figure 2. The data is also plotted.

- **`Figure3.R`**  
  Calculation of the data needed for Figure 3. The data is also plotted.

- **`Figure4.R`**  
  Calculation of the data needed for Figure 4. The data is also plotted.

- **`First_generation_probabilities.R`**  
  Probabilities for all cell types one generation after CRISPR was introduced and has taken full effect on the standing genetic variation (SGV).

- **`Numerical_approximation.R`**  
  Numerical calculation of the probability for extinction of a single cell lineage. De novo mutations are not included.

- **`Standing_genetic_variation.Rmd`**  
  Calculation of the SGV with Python (integrated in R). Must always be run before other files are executed.

- **`Supplemental_Figure8.R`**  
  Calculation of the data needed for Figure 8. The data is also plotted.

- **`Total_extinction_probability.R`**  
  Calculation of the extinction probability for the entire bacterial population.
