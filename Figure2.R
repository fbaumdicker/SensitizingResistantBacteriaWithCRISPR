

library(forcats)
library(ggplot2)
library(cowplot)


# Initiate parameters:

v <- 10^-8
s <- 0.1
SGV_Type <- "linear"
observed_copy_number_SGV <- 10
type_plasmid_duplication <- "random"
s_max <- 0
s_0 <- 0
u <- 0
dominance_function <- "linear"

# Set directory

main_directory <- dirname(rstudioapi::getSourceEditorContext()$path)

source(paste(main_directory, "/Total_extinction_population.R", sep = ""))

dir.create(gsub(" ", "", paste(main_directory, "/Images")), showWarnings = FALSE)

directory <- gsub(" ", "", paste(main_directory, "/Images/", as.character(Sys.Date())))

dir.create(directory, showWarnings = FALSE)

directory <- gsub(" ", "", paste(directory, "/"))


# Random replication probabilities:

calculate_probabilities_random <- function(number_of_pAMRmut_copies, observed_copy_number_SGV) { # RANDOM REPLICATION
  
  number_of_pCRISPR_copies <- 1
  
  if(number_of_pAMRmut_copies == observed_copy_number_SGV) {
    
    probabilities <- c(rep(0,observed_copy_number_SGV-1),1)
  
  } else {
    
    # Result vector to store probabilities
    probabilities <- numeric(observed_copy_number_SGV)
    
    # Calculate probabilities for each possible number of pAMRmut copies at the end
    for (k in 1:observed_copy_number_SGV) {
      if(k < number_of_pAMRmut_copies){
        probabilities[k] = 0
      }else{
        probabilities[k] <- ddirmnom(x = c(k-number_of_pAMRmut_copies,observed_copy_number_SGV-1-k), size = observed_copy_number_SGV-number_of_pAMRmut_copies-1, alpha = c(number_of_pAMRmut_copies, number_of_pCRISPR_copies), log = FALSE)
      }
    }
  }
    
  return(probabilities)
    
}

# Regular replication probabilities

calculate_probabilities_regular <- function(pAMRmut_copies, observed_copy_number_SGV) {  
  
  total_number_of_plasmids <- pAMRmut_copies + 1
  probabilities <- rep(0, observed_copy_number_SGV)
  
  while (total_number_of_plasmids < observed_copy_number_SGV/2) {
    pAMRmut_copies <- pAMRmut_copies * 2
    total_number_of_plasmids <- total_number_of_plasmids * 2
  }
  
  probabilities <- numeric(observed_copy_number_SGV)
  
  for (l in 0:(observed_copy_number_SGV-total_number_of_plasmids)) {
    m <- (observed_copy_number_SGV-total_number_of_plasmids) - l
    if (m > total_number_of_plasmids - pAMRmut_copies) next # Skip impossible cases
    
    probabilities[l + pAMRmut_copies] <- choose(pAMRmut_copies, l) * choose(total_number_of_plasmids - pAMRmut_copies, m) / choose(total_number_of_plasmids, (observed_copy_number_SGV-total_number_of_plasmids))
  }
  
  if (pAMRmut_copies == observed_copy_number_SGV) probabilities <- c(rep(0,observed_copy_number_SGV-1),1)
  
  return(probabilities)
}

SGV_data <- data.frame(copy_number = integer(), type = character(), frequency = double())

for (n in observed_copy_number_SGV) {
  
  SGV <- celltypefrequencies_R(v,n,s,SGV_Type)[1:n]
  # Compatible silencing:
  
  SGV_compatible_silencing <- SGV
  
  # Compatible cleaving:
  
  SGV_compatible_cleaving <- c(rep(0,n-1),sum(SGV_compatible_silencing[1:n]), SGV_compatible_silencing[n+1])
  
  # Incompatible silencing:
  
  SGV_incompatible_silencing <- SGV
  
  # Incompatible cleaving:
  
  SGV_incompatible_cleaving <- rep(0,n)
  
  for(i in 1:n){
    
    if (type_plasmid_duplication == "regular") {
      SGV_incompatible_cleaving <- SGV_incompatible_cleaving + calculate_probabilities_regular(i, observed_copy_number_SGV)*SGV[i]
    } else {
      SGV_incompatible_cleaving <- SGV_incompatible_cleaving + calculate_probabilities_random(i, observed_copy_number_SGV)*SGV[i]
    print(SGV_incompatible_cleaving)
      }
    
  }
  
  
  # Sum up the data
  
  for (i in 1:observed_copy_number_SGV) {
    SGV_data <- rbind(SGV_data, c(i, "SGV before pCRISPR introduction", SGV_compatible_silencing[i]))
  }
  
  for (i in 1:observed_copy_number_SGV) {
    SGV_data <- rbind(SGV_data, c(i, "Compatible cleaving", SGV_compatible_cleaving[i]))
  }
  
  for (i in 1:observed_copy_number_SGV) {
    SGV_data <- rbind(SGV_data, c(i, "Incompatible cleaving", SGV_incompatible_cleaving[i]))
  }
  
  for (i in 1:observed_copy_number_SGV) {
    SGV_data <- rbind(SGV_data, c(i, "Compatible silencing", SGV_compatible_silencing[i]))
  }
  
  for (i in 1:observed_copy_number_SGV) {
    SGV_data <- rbind(SGV_data, c(i, "Incompatible silencing", SGV_incompatible_silencing[i]))
  }
  
  
  colnames(SGV_data) <- c("copy_number", "type", "value")
  
  SGV_data$value <- as.numeric(SGV_data$value) 
  SGV_data$copy_number <- as.numeric(SGV_data$copy_number)
}


## Start the plot

SGV_data$type <- forcats::fct_relevel(SGV_data$type, "Compatible cleaving", "Incompatible cleaving", "Compatible silencing", "Incompatible silencing", "SGV before pCRISPR introduction")

plotA <- 
  ggplot(SGV_data, aes(x = copy_number, y = value, fill = factor(type))) +
  geom_bar(stat = "identity", position = "dodge") +
  suppressWarnings(scale_x_discrete(limits = c(1, 1:5*2))) +
  labs(title = "Sensitizing effect on the standing genetic variation", y = "Cell type frequency", x = "Number of mutated AMR plasmids") +
  scale_fill_manual(values = c("grey55", "darkorchid4", "cyan4", "darkorchid1", "cyan3"), limits = c("SGV before pCRISPR introduction", "Incompatible silencing", "Incompatible cleaving", "Compatible silencing", "Compatible cleaving")) +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black", size = 2),
        axis.text = element_text(size = 30, colour = "black"),
        axis.title = element_text(size = 30),
        title = element_text(size = 25),
        #legend.position = c(0.4,0.7),
        legend.position = "none",
        legend.text = element_text(size = 30), 
        axis.title.y = element_text(family = "serif", colour = "black", margin = margin(0, 15, 0, 0)),
        axis.title.x = element_text(family = "serif", margin = margin(15, 0, 0, 0)),
        text = element_text(size = 16, family = "serif"),
        plot.title = element_text(size = 35, face = "bold", hjust = 0.5, vjust = 4, margin = margin(30, 0, 0, 0)),
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"),
        legend.spacing.y = unit(0, "cm"),
        legend.key = element_rect(fill = NA),
        axis.ticks = element_line(size = 1.5),
        axis.ticks.length = unit(-0.2, "cm"),
        axis.text.x = element_text(margin = margin(10, 0, 0, 0)),  
        axis.text.y = element_text(margin = margin(0, 10, 0, 0)),
        axis.text.y.right = element_text(margin = margin(0, 0, 0, 10)))

 file <- paste(paste(directory, "SGV", ".png", sep = ""), sep = "")

ggsave(file = file,  units="px", width=720*5, height=480*5, dpi = 60*5, bg = "white")

#### LONGTERM SGV:

Longterm_SGV_data <- data.frame(cell_type = character(), type = character(), frequency = double())

for (cell_type in c("Compatible cleaving", "Compatible silencing", "Incompatible cleaving", "Incompatible silencing")) {
  Longterm_SGV_data <- rbind(Longterm_SGV_data, c("Homozygote mutant", cell_type, sum(subset(SGV_data, type == cell_type)$copy_number*subset(SGV_data, type == cell_type)$value)))
}

for (cell_type in c("Incompatible silencing", "Incompatible cleaving", "Compatible silencing", "Compatible cleaving")) {
  if (cell_type == "Incompatible cleaving") {
    Longterm_SGV_data <- rbind(Longterm_SGV_data, c("Homozygote CRISPR", cell_type, sum((observed_copy_number_SGV-subset(SGV_data, type == cell_type)$copy_number)*subset(SGV_data, type == cell_type)$value)))
  } else if (cell_type == "Incompatible silencing") {
    #Longterm_SGV_data <- rbind(Longterm_SGV_data, c("Homozygote CRISPR", cell_type, sum(1:observed_copy_number_SGV * rowSums(SGV_incompatible_silencing_matrix[1:observed_copy_number_SGV,1:observed_copy_number_SGV]))))
    Longterm_SGV_data <- rbind(Longterm_SGV_data, c("Homozygote CRISPR", cell_type, sum(SGV_compatible_silencing[1:observed_copy_number_SGV])*observed_copy_number_SGV/(observed_copy_number_SGV + 1)))
  } else {
    Longterm_SGV_data <- rbind(Longterm_SGV_data, c("Homozygote CRISPR", cell_type, 0))
  }
}

colnames(Longterm_SGV_data) <- c("cell_type", "type", "value")

Longterm_SGV_data$value <- as.numeric(Longterm_SGV_data$value)/observed_copy_number_SGV

Longterm_SGV_data$type <- forcats::fct_relevel(Longterm_SGV_data$type, 
                                               "Compatible cleaving", 
                                               "Incompatible cleaving", 
                                               "Compatible silencing", 
                                               "Incompatible silencing")

# Plot:

plotB <- 
  ggplot(Longterm_SGV_data, aes(x = cell_type, y = value, fill = factor(type))) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_x_discrete(limits = c("Homozygote mutant","Homozygote CRISPR")) +
  scale_fill_manual(values = c("darkorchid4", "cyan4", "darkorchid1", "cyan3"), limits = c("Incompatible silencing", "Incompatible cleaving", "Compatible silencing", "Compatible cleaving")) +
  labs(title = "Resistant population in the absence of antibiotics", y = "", x = "") +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black", size = 2),
        axis.text = element_text(size = 30, colour = "black"),
        axis.title = element_text(size = 30),
        title = element_text(size = 25),
        legend.position = "none",
        legend.text = element_text(size = 30), 
        axis.title.y = element_text(family = "serif", colour = "black", margin = margin(0, 15, 0, 0)),
        axis.title.x = element_text(family = "serif", margin = margin(15, 0, 0, 0)),
        text = element_text(size = 16, family = "serif"),
        plot.title = element_text(size = 35, face = "bold", hjust = 0.5, vjust = 4, margin = margin(30, 0, 0, 0)),
        plot.margin = margin(0, 40, 0, 0),
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"),
        legend.spacing.y = unit(0, "cm"),
        legend.key = element_rect(fill = NA),
        axis.ticks = element_line(size = 1.5),
        axis.ticks.length = unit(-0.2, "cm"),
        axis.text.x = element_text(margin = margin(10, 0, 0, 0)),  
        axis.text.y = element_text(margin = margin(0, 10, 0, 0)),
        axis.text.y.right = element_text(margin = margin(0, 0, 0, 10)))

file <- paste(paste(directory, "SGV_longterm", ".png", sep = ""), sep = "")

ggsave(file = file,  units="px", width=7200, height=4800, dpi = 600, bg = "white")



# Put the two plots together:

blank_space <- ggplot() + theme_void()

plot_grid(plotA, blank_space, plotB, 
          labels = c("A", "", "B"), 
          label_size = 35, 
          #label_fontface = "plain", 
          nrow = 1, 
          rel_widths = c(1, 0.1, 1))

ggsave(file = paste(paste(directory, "Fig2", ".png", sep = ""), sep = ""),  units="px", width=14400, height=4800, dpi = 600, bg = "white")

# The legend was added manually afterwards.

