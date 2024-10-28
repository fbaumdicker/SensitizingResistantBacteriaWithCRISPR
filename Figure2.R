

library(forcats)
library(ggplot2)
library(cowplot)
library(ggbreak)
library(ggpattern)
library(colorspace)

# Initiate parameters:

v <- 10^-8
s <- 0.1
SGV_Type <- "linear"
observed_copy_number_SGV <- 10
type_plasmid_duplication <- "random"
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
  
  SGV <- celltypefrequencies_R(v,n,s,SGV_Type)[c(n+1,1:n)]
  
  # Compatible silencing:
  
  SGV_compatible_silencing <- SGV
  
  # Compatible cleaving:
  
  SGV_compatible_cleaving <- c(SGV[1], rep(0,n-1),sum(SGV[2:(n+1)]))
  
  # Incompatible silencing:
  
  SGV_incompatible_silencing <- SGV
  
  # Incompatible cleaving:
  
  SGV_incompatible_cleaving <- c(SGV[1], rep(0,n))
  
  for(number_of_pAMRmut_copies in 1:n){
    
    if (type_plasmid_duplication == "regular") {
      SGV_incompatible_cleaving <- SGV_incompatible_cleaving + c(0, calculate_probabilities_regular(number_of_pAMRmut_copies, observed_copy_number_SGV)*SGV[number_of_pAMRmut_copies+1])
    } else {
      SGV_incompatible_cleaving <- SGV_incompatible_cleaving + c(0, calculate_probabilities_random(number_of_pAMRmut_copies, observed_copy_number_SGV)*SGV[number_of_pAMRmut_copies+1])
    }
    
  }
  
  
  # Sum up the data
  
  for (i in 1:(observed_copy_number_SGV+1)) {
    SGV_data <- rbind(SGV_data, c(i-1, "SGV before pCRISPR introduction", SGV_compatible_silencing[i]))
  }
  
  for (i in 1:(observed_copy_number_SGV+1)) {
    SGV_data <- rbind(SGV_data, c(i-1, "Compatible cleaving", SGV_compatible_cleaving[i]))
  }
  
  for (i in 1:(observed_copy_number_SGV+1)) {
    SGV_data <- rbind(SGV_data, c(i-1, "Incompatible cleaving", SGV_incompatible_cleaving[i]))
  }
  
  for (i in 1:(observed_copy_number_SGV+1)) {
    SGV_data <- rbind(SGV_data, c(i-1, "Compatible silencing", SGV_compatible_silencing[i]))
  }
  
  for (i in 1:(observed_copy_number_SGV+1)) {
    SGV_data <- rbind(SGV_data, c(i-1, "Incompatible silencing", SGV_incompatible_silencing[i]))
  }
  
  
  colnames(SGV_data) <- c("copy_number", "type", "value")
  
  SGV_data$value <- as.numeric(SGV_data$value) 
  SGV_data$copy_number <- as.numeric(SGV_data$copy_number)
}


## Start the plot

SGV_data$type <- forcats::fct_relevel(SGV_data$type, "Compatible cleaving", "Incompatible cleaving", "Compatible silencing", "Incompatible silencing", "SGV before pCRISPR introduction")

fill_colors <- c("SGV before pCRISPR introduction" = "grey55", 
                 "Incompatible silencing" = "darkorchid4", 
                 "Incompatible cleaving" = "cyan4", 
                 "Compatible silencing" = "darkorchid1", 
                 "Compatible cleaving" = "cyan3")

darkened_fill_colors <- sapply(fill_colors, function(color) darken(color, amount = 0.3))

plotA <-
  ggplot(SGV_data, aes(x = copy_number, y = value, fill = factor(type))) +
  geom_bar_pattern(stat = "identity", 
                   position = "dodge", 
                   pattern = ifelse(SGV_data$copy_number != 0, "stripe", "none"),
                   pattern_density = 0.3, 
                   pattern_angle = 45, 
                   pattern_spacing = 0.05,
                   pattern_fill = after_scale(darkened_fill_colors[as.character(SGV_data$type)]),
                   pattern_linetype = 1) +
  suppressWarnings(scale_x_discrete(limits = c(0, 1:5*2), labels = c("0 (WT)", 1:5*2))) +
  scale_y_continuous(breaks = c(c(1,3,5)*10^-7, 0.5,1), labels = c(c(1,3,5)*10^-7,0.5,1), limits = c(0,1.05)) +
  scale_y_break(c(6*10^-7,9*10^-7), scales = 0.6, ticklabels = NULL, expand = TRUE, space = 0.2) +
  labs(title = "Sensitizing effect on the standing genetic variation", y = "Cell type frequency", x = "Number of mutated AMR plasmids") +
  scale_fill_manual(values = fill_colors, 
                    limits = c("SGV before pCRISPR introduction", "Incompatible silencing", "Incompatible cleaving", "Compatible silencing", "Compatible cleaving")) +
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
        axis.title.y = element_text(angle = 90, family = "serif", colour = "black", margin = margin(0, 15, 0, 0)),
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
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank())

file <- paste(paste(directory, "SGV", ".png", sep = ""), sep = "")

ggsave(file = paste(paste(directory, "SGV", ".png", sep = ""), sep = ""),  units="px", width=7200, height=4800, dpi = 600, bg = "white")
 

#### LONGTERM SGV:

Longterm_SGV_data <- data.frame(cell_type = character(), type = character(), frequency = double())

for (cell_type in c("Incompatible silencing", "Incompatible cleaving", "Compatible silencing", "Compatible cleaving")) {
  
  if (cell_type == "Incompatible cleaving") {
    
    Longterm_SGV_data <- rbind(Longterm_SGV_data, c("Homozygote mutant due to mutants", cell_type, sum(subset(SGV_data, type == "SGV before pCRISPR introduction")$copy_number/(subset(SGV_data, type == "SGV before pCRISPR introduction")$copy_number+1)*subset(SGV_data, type == "SGV before pCRISPR introduction")$value)))
    Longterm_SGV_data <- rbind(Longterm_SGV_data, c("Homozygote CRISPR due to mutants", cell_type, sum(1/(subset(SGV_data, type == "SGV before pCRISPR introduction" & copy_number >= 1)$copy_number+1)*subset(SGV_data, type == "SGV before pCRISPR introduction" & copy_number >= 1)$value)))
    Longterm_SGV_data <- rbind(Longterm_SGV_data, c("Homozygote CRISPR due to WT cells", cell_type, subset(SGV_data, type == cell_type & copy_number == 0)$value ))
    
  } else if (cell_type == "Incompatible silencing") {
    
    Longterm_SGV_data <- rbind(Longterm_SGV_data, c("Homozygote mutant due to mutants", cell_type, 1/(observed_copy_number_SGV+1) * sum(subset(SGV_data, type == cell_type)$copy_number*subset(SGV_data, type == cell_type)$value)))
    Longterm_SGV_data <- rbind(Longterm_SGV_data, c("Homozygote CRISPR due to mutants", cell_type, 1/(observed_copy_number_SGV+1) *  sum(subset(SGV_data, type == cell_type & copy_number != 0)$value)))
    Longterm_SGV_data <- rbind(Longterm_SGV_data, c("Homozygote CRISPR due to WT cells", cell_type, 1/(observed_copy_number_SGV+1) *  sum(subset(SGV_data, type == cell_type & copy_number == 0)$value)))
    Longterm_SGV_data <- rbind(Longterm_SGV_data, c("Homozygote WT due to WT cells (non-silenced)", cell_type, observed_copy_number_SGV/(observed_copy_number_SGV+1) *  sum(subset(SGV_data, type == cell_type & copy_number == 0)$value)))
    Longterm_SGV_data <- rbind(Longterm_SGV_data, c("Homozygote WT due to mutants (non-silenced)", cell_type, 1/(observed_copy_number_SGV+1) * sum((observed_copy_number_SGV - subset(SGV_data, type == cell_type & copy_number != 0)$copy_number)*subset(SGV_data, type == cell_type & copy_number != 0)$value)))
    
  } else if (cell_type == "Compatible cleaving") {
    
    Longterm_SGV_data <- rbind(Longterm_SGV_data, c("Homozygote mutant due to mutants", cell_type, sum(subset(SGV_data, type == cell_type & copy_number >=1)$value)))
    Longterm_SGV_data <- rbind(Longterm_SGV_data, c("No plasmids (cleaved)", cell_type, subset(SGV_data, type == cell_type & copy_number == 0)$value ))
    
  } else if (cell_type == "Compatible silencing") {
    
    Longterm_SGV_data <- rbind(Longterm_SGV_data, c("Homozygote mutant due to mutants", cell_type, 1/observed_copy_number_SGV * sum(subset(SGV_data, type == cell_type)$copy_number*subset(SGV_data, type == cell_type)$value)))
    Longterm_SGV_data <- rbind(Longterm_SGV_data, c("Homozygote WT due to WT cells (silenced)", cell_type, subset(SGV_data, type == cell_type & copy_number == 0)$value ))
    Longterm_SGV_data <- rbind(Longterm_SGV_data, c("Homozygote WT due to mutants (silenced)", cell_type, 1/observed_copy_number_SGV *sum( (observed_copy_number_SGV - subset(SGV_data, type == cell_type & copy_number >= 1)$copy_number) *subset(SGV_data, type == cell_type & copy_number >= 1)$value)))
    
  }
}

colnames(Longterm_SGV_data) <- c("cell_type", "type", "value")

Longterm_SGV_data$value <- as.numeric(Longterm_SGV_data$value)

Longterm_SGV_data$type <- forcats::fct_relevel(Longterm_SGV_data$type, 
                                               "Compatible cleaving", 
                                               "Incompatible cleaving", 
                                               "Compatible silencing", 
                                               "Incompatible silencing")

Longterm_SGV_data_reduced <-  Longterm_SGV_data[Longterm_SGV_data$value > 0.01,]

for (i in 1:length(Longterm_SGV_data_reduced$cell_type)) {
  if (Longterm_SGV_data_reduced$cell_type[i] == "Homozygote CRISPR due to WT cells")
    Longterm_SGV_data_reduced$cell_type[i] <- "Homozygote CRISPR due to mutants"
  else if (Longterm_SGV_data_reduced$cell_type[i] == "Homozygote WT due to WT cells (non-silenced)")
    Longterm_SGV_data_reduced$cell_type[i] <- "Homozygote WT due to mutants (non-silenced)"
  else if (Longterm_SGV_data_reduced$cell_type[i] == "Homozygote WT due to WT cells (silenced)")
    Longterm_SGV_data_reduced$cell_type[i] <- "Homozygote WT due to mutants (silenced)"
}


plotB <-
  ggplot(Longterm_SGV_data, aes(x = cell_type, y = value, fill = factor(type))) +
  geom_bar(data = Longterm_SGV_data_reduced, stat = "identity", position = position_dodge2(preserve = "single"), width = 0.35) +
  geom_bar_pattern(data = Longterm_SGV_data[Longterm_SGV_data$value < 0.01,], 
                   aes(x = cell_type, y = value, fill = factor(type), 
                       pattern_fill = after_scale(darken(fill, 0.3))), 
                   stat = "identity", position = position_dodge2(preserve = "single"), 
                   pattern = "stripe", pattern_density = 0.3,  
                   pattern_spacing = 0.05, pattern_angle = 45, 
                   width = 0.7) +
  geom_vline(xintercept = 2.5, linetype = "dashed", color = "black", linewidth = 1) +
  annotate("text", x = 0.55, y = 0.91, label = "R", size = 14, fontface = "bold", color = "black", hjust = 0) + 
  annotate("text", x = 2.7, y = 0.91, label = "S", size = 14, fontface = "bold", hjust = 0.5) +
  scale_x_discrete(limits = c("Homozygote mutant due to mutants", "Homozygote WT due to mutants (non-silenced)", "Homozygote WT due to mutants (silenced)", "Homozygote CRISPR due to mutants", "No plasmids (cleaved)"), 
                   labels = c("hom. pAMRmut", "hom. pAMR\n(non-silenced)\n", "hom. pAMR \n(silenced)", "hom. pCRISPR", "no plasmid\n(cleaved pAMR)")) +
  scale_fill_manual(values = c("darkorchid4", "cyan4", "darkorchid1", "cyan3"), 
                    limits = c("Incompatible silencing", "Incompatible cleaving", "Compatible silencing", "Compatible cleaving")) +
  labs(title = "Resistant population in the absence of antibiotics", y = "", x = "") +
  scale_y_continuous(breaks = c(c(1,3,5)*10^-7, 0.5,1), labels = c(c(1,3,5)*10^-7,0.5,1), limits = c(0,1.05)) +
  scale_y_break(c(6*10^-7,9*10^-7), scales = 0.6, ticklabels = NULL, expand = TRUE, space = 0.2) +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black", linewidth = 2),
        axis.text = element_text(size = 30, colour = "black"),
        axis.title = element_text(size = 30),
        title = element_text(size = 25),
        legend.position = "none",
        legend.text = element_text(size = 30), 
        axis.title.y = element_text(family = "serif", colour = "black", margin = margin(0, 15, 0, 0)),
        axis.title.x = element_blank(),
        text = element_text(size = 16, family = "serif"),
        plot.title = element_text(size = 35, face = "bold", hjust = 0.5, vjust = 4, margin = margin(30, 0, 0, 0)),
        plot.margin = margin(0, 40, 0, 0),
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"),
        legend.spacing.y = unit(0, "cm"),
        legend.key = element_rect(fill = NA),
        axis.ticks = element_line(linewidth = 1.5),
        axis.ticks.length = unit(-0.2, "cm"),
        axis.ticks.y.right = element_blank(),
        axis.text.x = element_text(margin = margin(10, 0, 0, 0), size = 22), 
        axis.text.x.top =  element_blank(),  
        axis.text.y = element_text(margin = margin(0, 10, 0, 0)),
        axis.text.y.right = element_blank())


file <- paste(paste(directory, "SGV_longterm", ".png", sep = ""), sep = "")

ggsave(file = file,  units="px", width=7200, height=4800, dpi = 600, bg = "white")


# Afterwards the two plots were added together. The legend was later manually added. Some visuals were improved afterwards. 



