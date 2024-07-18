
library(ggplot2)

# Parameters:

v <- 10^-8
s <- 0.1
SGV_Type <- "linear"
observed_copy_numbers_supplement <- c(10,20,40)
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



# Calculation of the SGV for the different copy numbers:

SGV_supplement <- data.frame("Copy_Number" = sort(rep(observed_copy_numbers_supplement, max(observed_copy_numbers_supplement))), "Number_of_pAMRmut" = rep(1:max(observed_copy_numbers_supplement)/max(observed_copy_numbers_supplement), length(observed_copy_numbers_supplement)), "Frequency" = rep(0, max(observed_copy_numbers_supplement)*length(observed_copy_numbers_supplement)))
  
for (i in 1:length(observed_copy_numbers_supplement)) {
  
  n <- observed_copy_numbers_supplement[i]
  SGV_copy_number <- celltypefrequencies_R(v,n,s,SGV_Type)[1:n]
  
  for (k in 1:length(SGV_copy_number)) {
    SGV_supplement[SGV_supplement$Copy_Number == n & SGV_supplement$Number_of_pAMRmut == k/n, "Frequency"] <- SGV_copy_number[k]
  }
  
}

# Plot Supplemental Fig 8:
 
ggplot(SGV_supplement, aes(x = Number_of_pAMRmut, y = Frequency, fill = factor(Copy_Number))) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(title = "Standing genetic variation for increasing copy numbers", y = "Cell type frequency", x = "Fraction of mutated AMR plasmids") +
    scale_fill_manual(name = bquote("Copy numbers:"), values = c("blue4", "cornflowerblue", "aquamarine1"), labels = c("n = 10", "n = 20", "n = 40")) +
    theme(axis.line = element_line(color='black'),
          plot.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_rect(fill = "white", colour = "black", size = 2),
          axis.text = element_text(size = 30, colour = "black"),
          axis.title = element_text(size = 30),
          title = element_text(size = 25),
          legend.position = c(0.4,0.7),
          #legend.position = "none",
          legend.text = element_text(size = 30, margin = margin(l = 5, r = 10, t = 5, b = 5)), 
          axis.title.y = element_text(family = "serif", colour = "black", margin = margin(0, 15, 0, 0)),
          axis.title.x = element_text(family = "serif", margin = margin(15, 0, 0, 0)),
          text = element_text(size = 16, family = "serif"),
          plot.title = element_text(size = 35, face = "bold", hjust = 0.5, vjust = 4, margin = margin(30, 0, 0, 0)),
          legend.title = element_text(size = 30, family = "serif"),
          legend.background = element_blank(),
          legend.box.background = element_rect(colour = "black"),
          legend.margin = margin(20, 20, 20, 20),
          legend.spacing.y = unit(0.4, "cm"),
          legend.key = element_rect(fill = NA),
          axis.ticks = element_line(size = 1.5),
          axis.ticks.length = unit(-0.2, "cm"),
          axis.text.x = element_text(margin = margin(10, 0, 0, 0)),  
          axis.text.y = element_text(margin = margin(0, 10, 0, 0)),
          axis.text.y.right = element_text(margin = margin(0, 0, 0, 10)))
  
  file <- paste(paste(directory, "SGV_supplement", ".png", sep = ""), sep = "")
  
  ggsave(file = file,  units="px", width=9600, height=4800, dpi = 600, bg = "white")
  
