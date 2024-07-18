

library(reshape2)
library(ggplot2)
library(grid)
library(cowplot)


######### Parameters

# General: 

observed_copy_numbers <- 1:12
type_plasmid_duplication <- "regular"
u <- 0
number_of_mutations_at_start <- 1
used_cases <- 1:3
numerical_approximation_threshold <- 10^-10
SGV_Type <- "linear"
population_size <- 10^7
competition_ABR_plasmid <- 0
competition_CRISPR_plasmid <- 0

# Default (standard) cases

s_max_default <- 0.3
v_default <- 10^-8
s_default <- 0.1
s_0_default <- 0.2
dominance_function_default <- "linear"

# Varying parameters:

v_all <- c(10^-9, 10^(-7))
s_all <- c(0.01, 0.1, 0.9)
s_0_all <- c(0.05, 0.8)
s_max_all <- c(0.1, 1)
dominance_function_all <- c("dominant", "linear", "recessive")


###### Start:

main_directory <- dirname(rstudioapi::getSourceEditorContext()$path)

source(paste(main_directory, "/Total_extinction_population.R", sep = ""))

dir.create(gsub(" ", "", paste(main_directory, "/Images")), showWarnings = FALSE)

directory <- gsub(" ", "", paste(main_directory, "/Images/", as.character(Sys.Date())))

dir.create(directory, showWarnings = FALSE)

directory <- gsub(" ", "", paste(directory, "/"))


# List of variables which are to be observed. THIS HAS TO BE CHANGED TO PLOT PART A, B or C

variables <- c("s", "v", "s_max", "s_0", "dominance_function")
variables <- c("s") # Part A
# variables <- c("s", "v", "s_max", "s_0") # Part B
# variables <- c("dominance_function") # Part C

# Calculate point of intersection / the threshold copy number

find_intersection <- function(y1, y2) {
  for (i in 1:(length(y1)-1)) {
    if ((y1[i] - y2[i]) * (y1[i+1] - y2[i+1]) < 0 || (y2[i] - y1[i]) * (y2[i+1] - y1[i+1]) < 0) {
      x_intercept <- i - (y2[i] - y1[i]) / (y2[i+1] - y2[i] - y1[i+1] + y1[i])
      y_intercept <- y1[i] + (x_intercept-i) * (y1[i+1] - y1[i])
      return(c(x_intercept, y_intercept))
    }
  }
  return(NULL)  # No intersection point was found
}

# Initiate the list of all plots, which are put together in the end.

plotlist <- list()

# Go through all variables and compute every plot:

for (variable in variables) {
  
  if (variable == "s") {
    variable_all <- s_all
    SGV_Type <- "linear"
  } else if (variable == "v") {
    variable_all <- v_all
    SGV_Type <- "linear"
  } else if (variable == "s_0") {
    variable_all <- s_0_all
    SGV_Type <- "linear"
  } else if (variable == "s_max") {
    variable_all <- s_max_all
    SGV_Type <- "linear"
  } else if (variable == "dominance_function") {
    variable_all <- dominance_function_all
  }
  
  for (i in 1:length(variable_all)) {
    
    # Default values
    
    s_max <- s_max_default
    v <- v_default
    s <- s_default
    s_0 <- s_0_default
    dominance_function <- dominance_function_default
    
    # Variable: Overwrite certain value, which is observed
    
    if (variable == "s") {
      
      s <- s_all[i]
      variable_label <- bquote(paste(phantom(), "s" == -.(s_all[i])))
      
      if (i == 1) {
        title_image <- "High standing genetic variation"
      } else {
        title_image <- "Low standing genetic variation"
      }
      
    } else if (variable == "v") {
      v <- v_all[i]
      variable_label <- bquote(paste(phantom(), "N * v" == .(population_size*v_all[i])))
      
      if (i == 1) {
        title_image <- "Low mutational potential"
      } else {
        title_image <- "High mutational potential"
      }
      
    } else if (variable == "s_0") {
      s_0 <- s_0_all[i]
      variable_label <- bquote(paste(phantom(), lambda [min],  " = ", .(1-variable_all[i])))
      
      if (i == 1) {
        title_image <- "Weak antibiotic treatment"
      } else {
        title_image <- "Strong antibiotic treatment"
      }
      
    } else if (variable == "s_max") {
      s_max <- s_max_all[i]
      variable_label <- bquote(paste(phantom(), lambda [max],  " = ", .(1+variable_all[i])))
      if (i == 1) {
        title_image <- "Slow growth of mutant cells"
      } else {
        title_image <- "Fast growth of mutant cells"
      }
    
      } else if (variable == "dominance_function") {
      dominance_function <- dominance_function_all[i]
      SGV_Type <- "linear"
      variable_label <- ""
      if (i == 1) {
        title_image <- "Dominant fitness"
        y_axis_label <- "Total extinction probability"
      } else if (i == 2) {
        title_image <- "Linear fitness"
        y_axis_label <- ""
      } else if (i == 3) {
        title_image <- "Recessive fitness"
        y_axis_label <- ""
      }
    }
    
    if (variable == variables[[1]] && variable != "dominance_function") {
      y_axis_label <- "Total extinction probability"
    } else if (variable != "dominance_function"){
      y_axis_label <- ""
    }
 
    
    # Do the calculations of the rescue probabilities:
    
    main_result <- main(used_cases, observed_copy_numbers, u, dominance_function, number_of_mutations_at_start, type_plasmid_duplication, numerical_number_of_repeats, v, s_max, s_0, SGV_Type, population_size, competition_ABR_plasmid, competition_CRISPR_plasmid, numerical_approximation_threshold)
    rescue_probability_total_de_novo <- main_result[[4]]

    # Prepare the data for the plots:
    
    extinction_probability_total <- 1 - data.frame(rescue_probability_total_de_novo)
    colnames(extinction_probability_total) <- c("Compatible Type IV (silencing):\nWith de novo mutations", "Incompatible Type II (cleaving)", "Compatible Type II (cleaving)")
    extinction_probability_total <- cbind(Index = as.numeric(1:length(observed_copy_numbers)), extinction_probability_total)
    extinction_probability_total$Index <- as.integer(extinction_probability_total$Index)
    
    extinction_probability_total <- melt(extinction_probability_total, id.vars = 'Index', variable.name = 'series')
    extinction_probability_total$series <- factor(extinction_probability_total$series, levels = rev(levels( extinction_probability_total$series)))
    
    intersection_point <- find_intersection(1 - as.vector(rescue_probability_total_de_novo[[1]]), 1 - as.vector(rescue_probability_total_de_novo[[3]]))
    
    
    # Define the default case:
    if (i == 2 && length(variable_all) == 3 && variable != "dominance_function") {
      Is_default_case <- TRUE
    } else {
      Is_default_case <- FALSE
    }
    
    
    # Create a raster for the background color:
    
    steps <- 100
    custom_raster <- expand.grid(
      y = seq(-0.05, 1.05, length.out = steps),
      x = 1:12)
    custom_raster$color_value <- custom_raster$y
    custom_raster$color <- colorRampPalette(c("red","white"))(2*steps)[(steps + 1): (2*steps)]
    
    if (is.null(intersection_point)) {
      intersection_point_color <- custom_raster$color[steps + 1]
      intersection_point <- c(2,0)
      intersection_arrow_point <- c(-0.95, -1)
    } else {
      intersection_point_color <- "red2"
      intersection_arrow_point <- c(0.09, 0.14)
    }
    
    ymin <- max(-0.025, min(extinction_probability_total$value) - 0.05 * (max(extinction_probability_total$value) - min(extinction_probability_total$value)))
    ymax <- min(1.025,  max(extinction_probability_total$value) + 0.05 * (max(extinction_probability_total$value) - min(extinction_probability_total$value)))
    if (ymin < 0.15 & ymin >= 0) ymin <- 0
    if (ymax >= 0.85 & ymax <= 1) ymax <- 1
    
    custom_raster <- subset(custom_raster, y <= ymax & y >= ymin)
    
    # Plot:
    
    if (Is_default_case == FALSE) {
     plotlist[[length(plotlist)+1]] <- 
      ggdraw() + draw_plot(
        ggplot(extinction_probability_total, aes(x = Index, y = value)) + 
        geom_tile(data = custom_raster, aes(x = x, y = y, fill = color_value), alpha = 1) +
        coord_cartesian(expand = FALSE) +
        scale_fill_gradientn(colors = unique(custom_raster$color)) +
        geom_line(aes(colour = series), size = 2) +
        scale_colour_manual(values = c("cyan4", "cyan3", "darkorchid1", "darkorchid4")) +
        geom_point(aes(shape = series, color = series), size = 7) +
        geom_point(aes(x = intersection_point[1], y = intersection_point[2]), colour = intersection_point_color, pch = 4, size = 10, stroke = 3.5) +
        annotate("text", x = Inf, y = Inf, label = variable_label,
                 hjust = 1.25, vjust = 1.7, size = 16, family = "serif", color = "black") +
        scale_shape_manual(values = c(18, 16, 17, 15, 19)) +
        suppressWarnings(scale_x_discrete(limits = c(1, 1:6*2))) +
        scale_y_continuous(sec.axis = sec_axis(~.*1, name = ""), expand=c(0,0)) +
        labs(title = title_image, y = y_axis_label, x = expression(paste("Plasmid copy number ", italic("n")))) +
        theme(
          axis.line = element_line(color = 'black'),
          plot.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA, size = 2),
          panel.background = element_rect(colour = "black", size = 2),
          axis.text = element_text(size = 30, colour = "black"),
          axis.title = element_text(size = 30),
          title = element_text(size = 10),
          legend.position = "none",
          legend.text = element_text(size = 30),
          axis.title.y = element_text(family = "serif", colour = "black", margin = margin(0, 15, 0, 0)),
          axis.title.x = element_text(family = "serif", colour = "black", margin = margin(10, 0, 0, 0)),
          text = element_text(size = 20, family = "serif"),
          plot.title = element_text(size = 55, face = "bold", hjust = 0.5, vjust = 2, margin = margin(20, 0, 15, 0)),
          legend.title = element_blank(),
          legend.background = element_blank(),
          legend.box.background = element_rect(colour = "black"),
          legend.spacing.y = unit(0, "cm"),
          legend.key = element_rect(fill = NA),
          axis.ticks = element_line(size = 1.5),
          axis.ticks.length = unit(-0.2, "cm"),
          axis.text.x = element_text(margin = margin(10, 0, 0, 0)),
          axis.text.y = element_text(margin = margin(0, 10, 0, 0)),
          axis.text.y.right = element_text(margin = margin(0, 0, 0, 10))
        ) 
      ) + draw_plot(ggplot() + 
                      geom_segment(aes(x = intersection_point[1], y = intersection_arrow_point[1], xend = intersection_point[1], yend = intersection_arrow_point[2]),
                                   col = intersection_point_color,
                                   arrow = arrow(),
                                   size = 1.8,
                                   linejoin = "bevel",
                                   lineend = "round"
                      ) + theme_void() +
                      coord_cartesian(expand = FALSE, ylim = c(0,1), xlim = c(-1.14,13.93)), x = 0, y = 0)
      if (variable != "dominance_function") {
        if (i == 1) {
          file <- paste(paste(directory, "Observed_", sep = ""), "variable_", variable, paste("_","low", ".png", sep = ""), sep = "")
        } else if (i == 2 | i == 3) {
          file <- paste(paste(directory, "Observed_", sep = ""), "variable_", variable, paste("_","high", ".png", sep = ""), sep = "")
        } 
      } else if (variable == "dominance_function") {
          file <- paste(paste(directory, "Observed_", sep = ""), "variable_", variable, paste("_", variable_all[i] , ".png", sep = ""), sep = "")
      }
        
      ggsave(file = file,  units="px", width=7200, height=4800, dpi = 600, bg = "white")
      
    } else if (Is_default_case == TRUE) { # DEFAULT CASE
      
      title_image <- ""
      
      extinction_probability_total_typeII <- subset(extinction_probability_total, series %in% c("Incompatible Type II (cleaving)", "Compatible Type II (cleaving)"))
      extinction_probability_total_typeIV <- subset(extinction_probability_total, !(series %in% c("Incompatible Type II (cleaving)", "Compatible Type II (cleaving)")))
      colnames(extinction_probability_total_typeII)[2] <- "CRISPRTypeII"
      colnames(extinction_probability_total_typeIV)[2] <- "CRISPRTypeIV"
      
      # Plot the graph without legend:
      
      ggplot(extinction_probability_total, aes(x = Index, y = value)) + 
        geom_tile(data = custom_raster, aes(x = x, y = y, fill = color_value), alpha = 1) +
        coord_cartesian(expand = FALSE, ylim = c(0, 1)) +
        scale_shape_manual(name = "CRISPR Type IV (silencing):", values = c(18, 16, 17, 15, 19), labels = c("Incompatible with\nde novo mutations", "Incompatible without\nde novo mutations"), guide = guide_legend(order = 2)) +
        scale_fill_gradientn(name = "Extinction probability", colors = unique(custom_raster$color), limits = c(0,1), breaks = c(0,1), guide = guide_legend(order = 3)) +
        geom_line(aes(colour = series), size = 2) +
        scale_colour_manual(name = "CRISPR Type II (cleaving):", values = c("cyan4", "cyan3", "darkorchid1", "darkorchid4"), labels = c("Incompatible", "Compatible"), guide = guide_legend(order = 1)) +
        geom_point(aes(shape = series, color = series), size = 7) +
        geom_point(aes(x = intersection_point[1], y = intersection_point[2]), colour = intersection_point_color, pch = 4, size = 10, stroke = 3.5) +
        suppressWarnings(scale_x_discrete(limits = c(1, 1:6*2))) +
        labs(title = title_image, y = "Total extinction probability", x = expression(paste("Plasmid copy number ", italic("n")))) +
        theme(
          axis.line = element_line(color = 'black'),
          plot.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA, size = 2),
          panel.background = element_rect(colour = "black", size = 2),
          axis.text = element_text(size = 30, colour = "black"),
          axis.title = element_text(size = 30),
          title = element_text(size = 25),
          legend.position = "none",
          legend.text = element_text(size = 30, family = "serif"),
          axis.title.y = element_text(family = "serif", colour = "black", margin = margin(0, 80, 0, 0)),
          axis.title.x = element_text(family = "serif", margin = margin(15, 0, 0, 0)),
          text = element_text(size = 16, family = "serif"),
          plot.title = element_text(size = 60, face = "bold", hjust = 0.5, vjust = 2.5, margin = margin(30, 0, 0, 0)),
          legend.title = element_text(size = 30, family = "serif"),
          legend.background = element_blank(),
          legend.box.background = element_rect(colour = "black"),
          #legend.spacing.y = unit(3, "cm"),
          legend.key = element_rect(fill = NA),
          legend.margin = margin(-10, -10, 0, -10),
          axis.ticks = element_line(size = 1.5),
          axis.ticks.length = unit(-0.2, "cm"),
          axis.text.x = element_text(margin = margin(10, 0, 0, 0)),
          axis.text.y = element_text(margin = margin(0, 10, 0, 0)),
          axis.text.y.right = element_text(margin = margin(0, 0, 0, 10)),
          plot.margin = unit(c(0,10,0,2.5),"cm")
        ) +
        guides(
          fill = guide_colorbar(
            title = "Total extinction probability",
            title.position = "top",
            title.theme = element_text(family = "serif", size = 30),
            title.vjust = 1,
            barwidth = 1.5,
            barheight = 8,
            frame.colour = "black"
          )
        )
      
      
      file <- paste(paste(directory, "Default_case_plot", sep = ""), paste(".png", sep = ""), sep = "")
      
      
      ggsave(file = file,  units="px", width=9600, height=4800, dpi = 600, bg = "white")
      
      
      # Plot the legend: Part 1 (Needed only for the legend)
      
      ggplot(extinction_probability_total_typeII, aes(x = Index, y = value)) + 
        geom_tile(data = custom_raster, aes(x = x, y = y, fill = color_value), alpha = 1) +
        coord_cartesian(expand = FALSE, ylim = c(0, 1)) +
        scale_shape_manual(name = "pCleaving:", values = c(18, 16, 17, 15, 19), labels = c("Incompatible ", "Compatible"), guide = guide_legend(order = 1)) +
        scale_fill_gradientn(name = "Extinction probability", colors = unique(custom_raster$color), limits = c(0,1), breaks = c(0,1), guide = "none") +
        geom_line(aes(colour = CRISPRTypeII), size = 2) +
        scale_colour_manual(name = "pSilencing:", values = c("blue", "blue4", "darkorchid1", "darkorchid4"), labels = c("Compatible", ""), guide = guide_legend(order = 2), breaks = c("Compatible Type II (cleaving)", NA)) +
        geom_point(aes(shape = CRISPRTypeII, color = CRISPRTypeII), size = 7) +
        geom_point(aes(x = intersection_point[1], y = intersection_point[2]), colour = intersection_point_color, pch = 4, size = 10, stroke = 3.5) +
        suppressWarnings(scale_x_discrete(limits = c(1, 1:6*2))) +
        labs(title = title_image, y = y_axis_label, x = expression(paste("Plasmid copy number ", italic("n")))) +
        theme(
          axis.line = element_line(color = 'black'),
          plot.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA, size = 2),
          panel.background = element_rect(colour = "black", size = 2),
          axis.text = element_text(size = 30, colour = "black"),
          axis.title = element_text(size = 30),
          title = element_text(size = 25),
          legend.position = c(1.23, 0.537),
          legend.text = element_text(size = 30, family = "serif", margin=margin(t=7, b= 7)),
          axis.title.y = element_text(family = "serif", colour = "black", margin = margin(0, 80, 0, 0)),
          axis.title.x = element_text(family = "serif", margin = margin(15, 0, 0, 0)),
          text = element_text(size = 16, family = "serif"),
          plot.title = element_text(size = 50, face = "bold", hjust = 0.5, vjust = 2.5, margin = margin(30, 0, 0, 0)),
          legend.title = element_text(size = 30, family = "serif", margin=margin(t=0, b= 7)),
          legend.background = element_blank(),
          legend.box.background = element_rect(colour = "black"),
          legend.key = element_rect(fill = NA),
          legend.margin = margin(20, 10, 15, 20),
          legend.spacing.y = unit(0, "cm"),
          axis.ticks = element_line(size = 1.5),
          axis.ticks.length = unit(-0.2, "cm"),
          axis.text.x = element_text(margin = margin(10, 0, 0, 0)),
          axis.text.y = element_text(margin = margin(0, 10, 0, 0)),
          axis.text.y.right = element_text(margin = margin(0, 0, 0, 10)),
          plot.margin = unit(c(0,10,0,2.5),"cm")
        ) 
      
      file <- paste(paste(directory, "Default_case_legend_part_1", sep = ""), paste(".png", sep = ""), sep = "")
      
      ggsave(file = file,  units="px", width=9600, height=4800, dpi = 600, bg = "white")
      
      # Plot the legend: Part 2
      
      ggplot(extinction_probability_total_typeII, aes(x = Index, y = value)) + 
        geom_tile(data = custom_raster, aes(x = x, y = y, fill = color_value), alpha = 1) +
        coord_cartesian(expand = FALSE, ylim = c(0,1)) +
        scale_shape_manual(name = "Cleaving (CRISPR Type II):", values = c(18, 16, 17, 15, 19), labels = c("Incompatible", "Compatible"), guide = "none") +
        scale_fill_gradientn(name = "", colors = unique(custom_raster$color), limits = c(0,1), breaks = c(0,1), guide = guide_colourbar(direction = "horizontal",order = 3)) +
        geom_line(aes(colour = CRISPRTypeII), size = 2) +
        scale_colour_manual(name = "Silencing (CRISPR Type IV):", values = c("blue", "blue4", "darkorchid1", "darkorchid4"), labels = c("Compatible", ""), guide = "none") +
        geom_point(aes(shape = CRISPRTypeII, color = CRISPRTypeII), size = 7) +
        geom_point(aes(x = intersection_point[1], y = intersection_point[2]), colour = intersection_point_color, pch = 4, size = 10, stroke = 3.5) +
        suppressWarnings(scale_x_discrete(limits = c(1, 1:6*2))) +
        labs(title = title_image, y = y_axis_label, x = expression(paste("Plasmid copy number ", italic("n")))) +
        theme(
          axis.line = element_line(color = 'black'),
          plot.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA, size = 2),
          panel.background = element_rect(colour = "black", size = 2),
          axis.text = element_text(size = 30, colour = "black"),
          axis.title = element_text(size = 30),
          title = element_text(size = 25),
          legend.position = c(-0.15, 0.537),
          legend.text = element_text(size = 30, family = "serif"),
          axis.title.y = element_text(family = "serif", colour = "black", margin = margin(0, 80, 0, 0)),
          axis.title.x = element_text(family = "serif", margin = margin(15, 0, 0, 0)),
          text = element_text(size = 16, family = "serif"),
          plot.title = element_text(size = 50, face = "bold", hjust = 0.5, vjust = 2.5, margin = margin(30, 0, 0, 0)),
          legend.title = element_text(size = 30, family = "serif"),
          legend.background = element_blank(),
          legend.box.background = element_rect(colour = "white"),
          #legend.spacing.y = unit(3, "cm"),
          legend.key = element_rect(fill = NA),
          legend.margin = margin(-10, -10, 0, -10),
          axis.ticks = element_line(size = 1.5),
          axis.ticks.length = unit(-0.2, "cm"),
          axis.text.x = element_text(margin = margin(10, 0, 0, 0)),
          axis.text.y = element_text(margin = margin(0, 10, 0, 0)),
          axis.text.y.right = element_text(margin = margin(0, 0, 0, 10)),
          plot.margin = unit(c(0,10,0,2.5),"cm")
        ) +
        guides(
          fill = guide_colorbar(
            title = "",
            title.position = "top",
            title.theme = element_text(family = "serif", size = 30),
            barwidth = 1.5,
            barheight = 28.93,
            frame.colour = "black",
            label = FALSE,
          ) 
        )
      
      
      file <- paste(paste(directory, "Default_case_legend_part_2", sep = ""), paste(".png", sep = ""), sep = "")
      
      ggsave(file = file,  units="px", width=9600, height=4800, dpi = 600, bg = "white")
      
    } 
  }
}

## PART B:

blank_space <- ggplot() + theme_void()

separating_dash <- blank_space

modified_plotlist <- list(blank_space, blank_space, blank_space)

for (i in seq_along(plotlist)) {
  modified_plotlist <- c(modified_plotlist, plotlist[i])
  if (i %% 2 == 1 && i != length(plotlist)) {
    modified_plotlist <- c(modified_plotlist, list(blank_space))
  }
}

modified_plotlist_2 <- list()

for (i in seq_along(modified_plotlist)) {
  modified_plotlist_2 <- c(modified_plotlist_2, modified_plotlist[i])
  if (i %% 3 == 0 && i != length(modified_plotlist) && i != 3) {
    modified_plotlist_2 <- c(modified_plotlist_2, list(separating_dash), list(separating_dash), list(separating_dash))
  }
}

plot_grid(plotlist = modified_plotlist_2, 
          labels = c("", "B1", "", "B2", "", "B3", "", "B4"), 
          label_size = 60,
          label_fontfamily = "serif",
          label_x = -0.08,
          nrow = 3,
          ncol = 8,
          byrow = FALSE,
          rel_widths = c(0.05, 1, 0.05, 1, 0.05, 1, 0.05, 1),
          rel_heights = c(1, 0.1 ,1))

ggsave(file = paste(paste(directory, "Fig4B", ".png", sep = ""), sep = ""),  units="px", width=2880*2.1, height=960*2.1, dpi = 60*2.1, bg = "white", limitsize = FALSE)

## PART C:

modified_plotlist_C <- list()

for (i in seq_along(plotlist)) {
  modified_plotlist_C <- c(modified_plotlist_C, list(blank_space), plotlist[i])
  }

plot_grid(plotlist = modified_plotlist_C, 
          labels = c("C1", "", "C2", "", "C3"), 
          label_size = 55,
          label_fontfamily = "serif",
          label_x = -0.05,
          nrow = 1,
          ncol = 6,
          rel_widths = c(0.03, 1, 0.03, 1, 0.03, 1)
          )

ggsave(file = paste(paste(directory, "Fig4C", ".png", sep = ""), sep = ""),  units="px", width=960*3*0.9, height=480*1/0.9, dpi = 60*1/0.9, bg = "white")


