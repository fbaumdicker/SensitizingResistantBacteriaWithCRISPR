

library(reshape2)
library(ggplot2)
library(cowplot)


######### Parameters

# General: 

observed_copy_numbers <- 1:12
type_plasmid_duplication <- "random"
dominance_function <- "linear"
SGV_Type <- "linear"
u <- 0
numerical_approximation_threshold <- 10^-10
used_cases <- 1:6
competition_ABR_plasmid <- 0
competition_CRISPR_plasmid <- 0

# Default cases

lambda_max <- 0.05
v <- 10^-8
s <- 0.1
population_size <- 10^7
lambda_min <- 0.05


###### Start:

# Set directory

main_directory <- dirname(rstudioapi::getSourceEditorContext()$path)

source(paste(main_directory, "/Total_extinction_population.R", sep = ""))

dir.create(gsub(" ", "", paste(main_directory, "/Images")), showWarnings = FALSE)

directory <- gsub(" ", "", paste(main_directory, "/Images/", as.character(Sys.Date())))

dir.create(directory, showWarnings = FALSE)

directory <- gsub(" ", "", paste(directory, "/"))



###### Plot 1: Wild-type extinction probability

number_of_mutations_at_start <- 0

main_result_wild_type <- main(used_cases, observed_copy_numbers, u, dominance_function, number_of_mutations_at_start, type_plasmid_duplication, numerical_number_of_repeats, v, lambda_max, lambda_min, SGV_Type, population_size, competition_ABR_plasmid, competition_CRISPR_plasmid, numerical_approximation_threshold)

extermination_probability <- main_result_wild_type[[1]]

# Prepare the data:

extermination_probability_df <- 1-data.frame(extermination_probability)
colnames(extermination_probability_df) <- c("Compatible Silencing", "Compatible Cleaving", "Incompatible Cleaving", "Incompatible Silencing", "Compatible No Effect", "Incompatible No Effect")
extermination_probability_df <- cbind(Index = as.numeric(1:length(observed_copy_numbers)), extermination_probability_df)
extermination_probability_df$Index <- as.integer(extermination_probability_df$Index)

extermination_probability_df <- melt(extermination_probability_df ,  id.vars = 'Index', variable.name = 'series')
extermination_probability_df$series <- factor(extermination_probability_df$series, levels = rev(levels( extermination_probability_df$series)))

# Plot: Show only every third point for points with certain extermination:

extermination_probability_df <- subset(extermination_probability_df, series %in% c("Incompatible Silencing", "Incompatible No Effect", "Compatible No Effect") | (series == "Compatible Silencing" & Index%%3 == 0) | (series == "Incompatible Cleaving" & Index%%3 == 1) | (series == "Compatible Cleaving" & Index%%3 == 2))

title_image <- "Single wild-type cell"

plot_A <- 
  ggplot(extermination_probability_df, aes(x = Index, y = value)) + 
  geom_line(aes(colour = series), size = 2) +
  geom_segment(aes(x = 1, y = 1, xend = 12, yend = 1), color = "cyan4", size = 2) +
  scale_colour_manual(values =  c("grey30","grey55", "darkorchid4", "cyan4", "cyan3","darkorchid1")) +
  geom_point(aes(shape = series, color = series), size = 7) +
  scale_shape_manual(values = c(18, 20, 15, 16, 19, 17)) +
  suppressWarnings(scale_x_discrete(limits = c(1, 1:6*2))) +
  labs(title = title_image, y = "Single-cell extinction probability", x = expression(paste("Plasmid copy number ", italic("n")))) +
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



file <- paste(paste(directory, "Single_wild_type_cell_plot", sep = ""), paste(".png", sep = ""), sep = "")


ggsave(file = file,  units="px", width = 7200, height=4800, dpi = 600, bg = "white")

# Plot the legend:

extermination_probability_df_legend <- subset(extermination_probability_df,  series %in% c( "Compatible No Effect", "Compatible Silencing"))

ggplot(extermination_probability_df_legend, aes(x = Index, y = value)) + 
  geom_line(aes(colour = series), size = 2) +
  scale_colour_manual(name = bquote("pCleaving:"), values =  c("grey30","grey55", "darkorchid4", "cyan4", "cyan3","darkorchid1"), guide = guide_legend(order = 1), labels = c("Compatible", "Incompatible")) +
  geom_point(aes(shape = series, color = series, alpha = series), size = 7) +
  scale_shape_manual(name = bquote("pSilencing:"), values = c(18, 20, 15, 16, 19, 17), guide = guide_legend(order = 2), labels = c("Compatible", "Incompatible")) +
  scale_alpha_manual(values = c(1.5, 1.5), name = bquote(paste("Non", "-", "CRISPR plasmids:")), guide = guide_legend(order = 3), labels = c("Incompatible", paste("Compatible/No effect"))) +
  suppressWarnings(scale_x_discrete(limits = c(1, 1:6*2))) +
  labs(title = title_image, y = "", x = expression(paste("Plasmid copy number ", italic("n")))) +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black", size = 2),
        axis.text = element_text(size = 30, colour = "black"),
        axis.title = element_text(size = 30),
        title = element_text(size = 25),
        legend.position = c(0.7, 0.5),
        legend.text = element_text(size = 25), 
        axis.title.y = element_text(family = "serif", colour = "black", margin = margin(0, 15, 0, 0)),
        axis.title.x = element_text(family = "serif", margin = margin(15, 0, 0, 0)),
        text = element_text(size = 16, family = "serif"),
        plot.title = element_text(size = 35, face = "bold", hjust = 0.5, vjust = 4, margin = margin(30, 0, 0, 0)),
        legend.title = element_text(size = 25, family = "serif"),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"),
        legend.spacing.y = unit(0, "cm"),
        legend.key = element_rect(fill = NA),
        axis.ticks = element_line(size = 1.5),
        axis.ticks.length = unit(-0.2, "cm"),
        axis.text.x = element_text(margin = margin(10, 0, 0, 0)),  
        axis.text.y = element_text(margin = margin(0, 10, 0, 0)),
        axis.text.y.right = element_text(margin = margin(0, 0, 0, 10))
  ) + guides(alpha = guide_legend(override.aes = list()))


file <- paste(paste(directory, "Single_wild_type_cell_legend", sep = ""), paste(".png", sep = ""), sep = "")


ggsave(file = file,  units="px", width = 7200, height=4800, dpi = 600, bg = "white")



###### Plot 2: Single mutated cell

number_of_mutations_at_start <- 1

main_result_single_cell <- main(used_cases, observed_copy_numbers, u, dominance_function, number_of_mutations_at_start, type_plasmid_duplication, numerical_number_of_repeats, v, lambda_max, lambda_min, SGV_Type, population_size, competition_ABR_plasmid, competition_CRISPR_plasmid, numerical_approximation_threshold)

extermination_probability <- main_result_single_cell[[1]]

# Prepare the data:

extermination_probability_df <- 1-data.frame(extermination_probability)
colnames(extermination_probability_df) <- c("Compatible Silencing", "Compatible Cleaving", "Incompatible Cleaving", "Incompatible Silencing", "Compatible No Effect", "Incompatible No Effect")
extermination_probability_df <- cbind(Index = as.numeric(1:length(observed_copy_numbers)), extermination_probability_df)
extermination_probability_df$Index <- as.integer(extermination_probability_df$Index)

extermination_probability_df <- melt(extermination_probability_df ,  id.vars = 'Index', variable.name = 'series')
extermination_probability_df$series <- factor(extermination_probability_df$series, levels = rev(levels( extermination_probability_df$series)))

extermination_probability_df_one_mutation <- subset(extermination_probability_df, series %in% c("Incompatible Silencing", "Incompatible No Effect", "Compatible Silencing", "Incompatible Cleaving") | (series == "Compatible No Effect" & Index%%2 == 0) | (series == "Compatible Cleaving" & Index%%2 == 1))

title_image <- "Cell with one mutated AMR plasmid"

plot_B <- 
  ggplot(extermination_probability_df_one_mutation, aes(x = Index, y = value)) + 
  geom_line(aes(colour = series), size = 2) +
  geom_segment(aes(x = 1, y = min(value), xend = 12, yend = min(value)), color = "grey55", size = 2) +
  scale_colour_manual(values =  c("grey30","grey55", "darkorchid4", "cyan4", "cyan3","darkorchid1")) +
  geom_point(aes(shape = series, color = series), size = 7) +
  scale_shape_manual(values = c(18, 20, 15, 16, 19, 17)) +
  suppressWarnings(scale_x_discrete(limits = c(1, 1:6*2))) +
  labs(title = title_image, y = "", x = expression(paste("Plasmid copy number ", italic("n")))) +
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
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"),
        legend.spacing.y = unit(0, "cm"),
        legend.key = element_rect(fill = NA),
        axis.ticks = element_line(size = 1.5),
        axis.ticks.length = unit(-0.2, "cm"),
        axis.text.x = element_text(margin = margin(10, 0, 0, 0)),  
        axis.text.y = element_text(margin = margin(0, 10, 0, 0)),
        axis.text.y.right = element_text(margin = margin(0, 0, 0, 10))) +
  coord_cartesian(ylim = c(min(extermination_probability_df_legend$value),1))



file <- paste(paste(directory, "Cell_with_one_mutated_plasmid", sep = ""), paste(".png", sep = ""), sep = "")


ggsave(file = file,  units="px", width=7200, height=4800, dpi = 600, bg = "white")


# Put the two plots together. The legend was later manually added. Some visuals were improved afterwards. 

blank_space <- ggplot() + theme_void()

plot_grid(plot_A, blank_space, plot_B, 
          labels = c("A", "", "B"), 
          label_size = 35, 
          nrow = 1, 
          rel_widths = c(1, 0.1, 1))

ggsave(file = paste(paste(directory, "Fig3", ".png", sep = ""), sep = ""),  units="px", width=14400, height=4800, dpi = 600, bg = "white")






