

library(reshape2)
library(ggplot2)
library(grid)
library(cowplot)


n <- 12
copy_numbers <- 0:n
lambda_min <- 0.8
lambda_max_all <- c(1.1, 1.3, 2)

# Initialize:

plotlist <- list()

# Calculate birth-death ratios:

birth_death_ratio_dominant_lambda_1 <- lambda_min + (lambda_max_all[1]-lambda_min)*c(0, rep(1, n))
birth_death_ratio_linear_lambda_1 <- lambda_min + (lambda_max_all[1]-lambda_min)*0:n/n
birth_death_ratio_recessive_lambda_1 <- lambda_min + (lambda_max_all[1]-lambda_min)*c(rep(0, n), 1)

birth_death_ratio_dominant_lambda_2 <- lambda_min + (lambda_max_all[2]-lambda_min)*c(0, rep(1, n))
birth_death_ratio_linear_lambda_2 <- lambda_min + (lambda_max_all[2]-lambda_min)*0:n/n
birth_death_ratio_recessive_lambda_2 <- lambda_min + (lambda_max_all[2]-lambda_min)*c(rep(0, n), 1)

birth_death_ratio_dominant_lambda_3 <- lambda_min + (lambda_max_all[3]-lambda_min)*c(0, rep(1, n))
birth_death_ratio_linear_lambda_3 <- lambda_min + (lambda_max_all[3]-lambda_min)*0:n/n
birth_death_ratio_recessive_lambda_3 <- lambda_min + (lambda_max_all[3]-lambda_min)*c(rep(0, n), 1)


birth_death_ratio_dominant_df <- data.frame(rbind(birth_death_ratio_dominant_lambda_1,birth_death_ratio_dominant_lambda_2,birth_death_ratio_dominant_lambda_3))
rownames(birth_death_ratio_dominant_df) <- c("lambda_1", "lambda_2", "lambda_3")
colnames(birth_death_ratio_dominant_df) <- 0:n
birth_death_ratio_dominant_df <- melt(as.matrix(birth_death_ratio_dominant_df),varnames = c("series", "Index"),value.name = "value")
birth_death_ratio_dominant_df$Index <- as.numeric(as.character(birth_death_ratio_dominant_df$Index)) 

birth_death_ratio_linear_df <- data.frame(rbind(birth_death_ratio_linear_lambda_1,birth_death_ratio_linear_lambda_2,birth_death_ratio_linear_lambda_3))
rownames(birth_death_ratio_linear_df) <- c("lambda_1", "lambda_2", "lambda_3")
colnames(birth_death_ratio_linear_df) <- 0:n
birth_death_ratio_linear_df <- melt(as.matrix(birth_death_ratio_linear_df),varnames = c("series", "Index"),value.name = "value")
birth_death_ratio_linear_df$Index <- as.numeric(as.character(birth_death_ratio_linear_df$Index)) 

birth_death_ratio_recessive_df <- data.frame(rbind(birth_death_ratio_recessive_lambda_1,birth_death_ratio_recessive_lambda_2,birth_death_ratio_recessive_lambda_3))
rownames(birth_death_ratio_recessive_df) <- c("lambda_1", "lambda_2", "lambda_3")
colnames(birth_death_ratio_recessive_df) <- 0:n
birth_death_ratio_recessive_df <- melt(as.matrix(birth_death_ratio_recessive_df),varnames = c("series", "Index"),value.name = "value")
birth_death_ratio_recessive_df$Index <- as.numeric(as.character(birth_death_ratio_recessive_df$Index)) 


for (df in list(birth_death_ratio_recessive_df, birth_death_ratio_linear_df, birth_death_ratio_dominant_df)) {
  if (df$value[4] == lambda_max_all[1]) {
    #legend_position = c(0.05, 0.8)
    x_axis_title = expression(paste("Number of functional AMR plasmids ", italic(k)))
  } else {
    #legend_position = "none"
    x_axis_title = ""
  } 
    
  if (df$value[4] == lambda_min) {
     df <- rbind(
      subset(df, series == "lambda_1" & (Index %% 3 == 0 | Index == n)),
      subset(df, series == "lambda_2" & (Index %% 3 == 1 | Index == n)),
      subset(df, series == "lambda_3" & (Index %% 3 == 2 | Index == n))
    )
    legend_position = c(0.28, 0.85)
    title_first_column = bquote(bold("Birth-death ratio " * italic("R")))
  } else {

    legend_position = "none"
    title_first_column = ""
  }
    
  plotlist[[length(plotlist)+1]] <- 
    ggplot(df, aes(x = Index, y = value)) + 
      geom_point(aes(color = series), shape = 18, size = 5, stroke = 4) +
      scale_colour_manual(values = c("cyan3","darkorchid4", "grey30", "grey55", "darkorchid4", "cyan4", "cyan3", "darkorchid1"), 
                          labels = c(expression(paste(" ", italic(R)[max], " = 1.1")), expression(paste(" ", italic(R)[max], " = 1.3")), expression(paste(" ", italic(R)[max], " = 2")))) +
                          suppressWarnings(scale_x_discrete(limits = c(0, 1:6*2))) +
      labs(title = title_first_column,
           y = bquote(italic("R")[k] * " = " * lambda[k] * " / " * mu[k]),
           x = x_axis_title) +
      theme(axis.line = element_line(color = 'black'),
            plot.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_rect(fill = "white", colour = "black", size = 2),
            axis.text = element_text(size = 30, colour = "black"),
            axis.title = element_text(size = 30),
            title = element_text(size = 25),
            legend.position = legend_position,
            legend.justification = c(0, 1),
            legend.text = element_text(size = 30, margin = margin(b = 3)), 
            axis.title.y = element_text(family = "serif", colour = "black", margin = margin(0, 15, 0, 0)),
            axis.title.x = element_text(family = "serif", margin = margin(15, 0, 0, 0)),
            text = element_text(size = 16, family = "serif"),
            plot.title = element_text(size = 35, face = "bold", hjust = 0.5, vjust = 4, margin = margin(30, 0, 0, 0)),
            legend.title = element_blank(),
            legend.background = element_blank(),
            legend.box.background = element_rect(colour = "black"),
            legend.key = element_rect(fill = NA, colour = NA),
            legend.key.width = unit(1.1, "cm"),
            legend.margin = margin(10, 15, 10, 15),
            axis.ticks = element_line(size = 1.5),
            axis.ticks.length = unit(-0.2, "cm"),
            axis.text.x = element_text(margin = margin(10, 0, 0, 0)),  
            axis.text.y = element_text(margin = margin(0, 10, 0, 0)),
            axis.text.y.right = element_text(margin = margin(0, 0, 0, 10))) +
      coord_cartesian(ylim = c(min(birth_death_ratio_recessive_df$value)-0.05, 
                               max(birth_death_ratio_recessive_df$value)+0.05)) +
      guides(
        color = guide_legend(
          override.aes = list(size = 5),
          byrow = TRUE
        )
      )
}


#### Birth and death parameters:

# Bacteriostatic drug:

birth_death_parameters_bacteriostatic <- function(df) {
  df_mu <- df
  df_mu$series <- sub("lambda", "mu", df_mu$series)
  df_mu$value <- 1
  
  return(df)
}

birth_death_parameters_bacteriostatic_dominant <- birth_death_parameters_bacteriostatic(birth_death_ratio_dominant_df)
birth_death_parameters_bacteriostatic_linear <- birth_death_parameters_bacteriostatic(birth_death_ratio_linear_df)
birth_death_parameters_bacteriostatic_recessive <- birth_death_parameters_bacteriostatic(birth_death_ratio_recessive_df)

# Bactericidal drug:

birth_death_parameters_bactericidal <- function(df) {
  df_lambda <- df
  df_lambda$value <- 1
  df$series <- gsub("lambda", "mu", df$series)
  df$value <- 1 / df$value
  
  return(df)
}

birth_death_parameters_bactericidal_dominant <- birth_death_parameters_bactericidal(birth_death_ratio_dominant_df)
birth_death_parameters_bactericidal_linear <- birth_death_parameters_bactericidal(birth_death_ratio_linear_df)
birth_death_parameters_bactericidal_recessive <- birth_death_parameters_bactericidal(birth_death_ratio_recessive_df)


# Plots

for (df in list(birth_death_ratio_recessive_df, birth_death_ratio_linear_df, birth_death_ratio_dominant_df)) {
  if (df$value[4] == lambda_max_all[1]) {
    x_axis_title = expression(paste("Number of functional AMR plasmids ", italic("k")))
  } else {
    x_axis_title = ""
    legend_position = "none"
  }
  if (df$value[4] == lambda_min) {
    plot_title = "Bacteriostatic drug"
    legend_position = "none"
    text = bquote(mu[k] %==%  1)
    df <- rbind(
      subset(df, series == "lambda_1" & (Index %% 3 == 0 | Index == n)),
      subset(df, series == "lambda_2" & (Index %% 3 == 1 | Index == n)),
      subset(df, series == "lambda_3" & (Index %% 3 == 2 | Index == n))
    )
  } else {
    plot_title = ""
    text = ""
    legend_position = "none"
  }
  
  df_plot <- birth_death_parameters_bacteriostatic(df)
  plotlist[[length(plotlist)+1]] <- 
    ggplot(df_plot, aes(x = Index, y = value)) + 
      geom_point(aes(color = series, shape = series), size = 3, stroke = 4) +
      scale_shape_manual(values = c(
        "lambda_1" = 3, "lambda_2" = 3, "lambda_3" = 3,
        "mu_1" = 19, "mu_2" = 19, "mu_3" = 19
      ), labels = c(
        bquote(" Birth rates " * lambda[k] * " for " * lambda[max] ~ "=" ~ .(lambda_max_all[1])), bquote(" Birth rates " * lambda[k] * " for " * italic(R)[max] ~ "=" ~ .(lambda_max_all[2])), bquote(" Birth rates " * lambda[k] * " for " * italic(R)[max] ~ "=" ~ .(lambda_max_all[3])), expression(" Death rate"~mu[k] %==% 1)
      )) +
      geom_hline(yintercept = 1, linetype = "dashed", color = "black", size = 2, alpha = 0.7) +
      scale_colour_manual(values = c("cyan3","darkorchid4", "grey30", "darkred", "darkorchid4", "cyan4", "cyan3", "darkorchid1"), 
                          labels = c(expression(paste(" ", italic("R"), " = 1.1")), expression(paste(" ", italic("R"), " = 1.3")), expression(paste(" ", italic("R"), " = 2")))) +
      suppressWarnings(scale_x_discrete(limits = c(0, 1:6*2))) +
      labs(title = plot_title,
           y = bquote("Birth rates " * lambda[k]),
           x = x_axis_title) +
    annotate("text", 
             x = min(df_plot$Index) + 1.5,     
             y = min(df_plot$value) + 0.4,    
             label = text,        
             hjust = 0.5,            
             vjust = 0.5,                       
             size = 12,                         
             family = "serif", 
             fontface = "italic", 
             color = "black",
             alpha = 0.7) +
      theme(axis.line = element_line(color = 'black'),
            plot.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_rect(fill = "white", colour = "black", size = 2),
            axis.text = element_text(size = 30, colour = "black"),
            axis.title = element_text(size = 30),
            title = element_text(size = 25),
            legend.position = legend_position,
            legend.justification = c(0, 1),
            legend.text = element_text(size = 30, margin = margin(b = 3)), 
            axis.title.y = element_text(family = "serif", colour = "black", margin = margin(0, 15, 0, 0)),
            axis.title.x = element_text(family = "serif", margin = margin(15, 0, 0, 0)),
            text = element_text(size = 16, family = "serif"),
            plot.title = element_text(size = 35, face = "bold", hjust = 0.5, vjust = 4, margin = margin(30, 0, 0, 0)),
            legend.title = element_blank(),
            legend.background = element_blank(),
            legend.box.background = element_rect(colour = "black"),
            legend.key = element_rect(fill = NA, colour = NA),
            legend.key.width = unit(1.1, "cm"),
            legend.margin = margin(10, 15, 10, 15),
            axis.ticks = element_line(size = 1.5),
            axis.ticks.length = unit(-0.2, "cm"),
            axis.text.x = element_text(margin = margin(10, 0, 0, 0)),  
            axis.text.y = element_text(margin = margin(0, 10, 0, 0)),
            axis.text.y.right = element_text(margin = margin(0, 0, 0, 10))) +
      coord_cartesian(ylim = c(min(df_plot$value)-0.05, 
                               max(df_plot$value)+0.05)) +
      guides(
        color = "none"
      )
}  


for (df in list(birth_death_ratio_recessive_df, birth_death_ratio_linear_df, birth_death_ratio_dominant_df)) {
  if (df$value[4] == lambda_max_all[1]) {
    legend_position = "none"
    x_axis_title = expression(paste("Number of functional AMR plasmids ", italic(k)))
  } else {
    legend_position = "none"
    x_axis_title = ""
  }
  if (df$value[4] == lambda_min) {
    df <- rbind(
      subset(df, series == "lambda_1" & (Index %% 3 == 0 | Index == n)),
      subset(df, series == "lambda_2" & (Index %% 3 == 1 | Index == n)),
      subset(df, series == "lambda_3" & (Index %% 3 == 2 | Index == n))
    )
    plot_title = "Bactericidal drug"
    text = bquote(lambda[k] %==%  1)
  } else {
    plot_title = ""
    text = ""
  }

  
  df_plot <- birth_death_parameters_bactericidal(df)
  plotlist[[length(plotlist)+1]] <- 
    ggplot(df_plot, aes(x = Index, y = value)) + 
    geom_point(aes(color = series, shape = series), size = 3, stroke = 4) +
    scale_shape_manual(values = c(
      "lambda_1" = 3, "lambda_2" = 3, "lambda_3" = 3,
      "mu_1" = 19, "mu_2" = 19, "mu_3" = 19
    ), labels = c(
      bquote(" Death rates " * mu[k] * " for " * italic(R)[max] ~ "=" ~ .(lambda_max_all[1])), bquote(" Death rates " * mu[k] * " for " * italic(R)[max] ~ "=" ~ .(lambda_max_all[2])), bquote(" Death rates " * mu[k] * " for " * italic(R)[max] ~ "=" ~ .(lambda_max_all[3])))
    ) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "black", size = 2, alpha = 0.7) +
    scale_colour_manual(values = c( #"darkgreen", 
      "cyan3","darkorchid4", "grey30", "black", "darkorchid4", "cyan4", "cyan3", "darkorchid1"), 
                        labels = c(expression(paste(" ", italic("R"), " = 1.1")), expression(paste(" ", italic("R"), " = 1.3")), expression(paste(" ", italic("R"), " = 2")))) +
    suppressWarnings(scale_x_discrete(limits = c(0, 1:6*2))) +
    labs(title = plot_title,
         y = bquote("Death rates " * mu[k]),
         x = x_axis_title) +
    annotate("text", 
               x = min(df_plot$Index) + 1.5,     
               y = max(df_plot$value) - 0.4,    
               label = text,        
               hjust = 0.5,                       
               vjust = 0.5,                      
               size = 12,                       
               family = "serif", 
               fontface = "italic", 
               color = "black",
               alpha = 0.7) +
    theme(axis.line = element_line(color = 'black'),
          plot.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_rect(fill = "white", colour = "black", size = 2),
          axis.text = element_text(size = 30, colour = "black"),
          axis.title = element_text(size = 30),
          title = element_text(size = 25),
          legend.position = legend_position,
          legend.justification = c(0, 1),
          legend.text = element_text(size = 30, margin = margin(b = 3)), 
          axis.title.y = element_text(family = "serif", colour = "black", margin = margin(0, 15, 0, 0)),
          axis.title.x = element_text(family = "serif", margin = margin(15, 0, 0, 0)),
          text = element_text(size = 16, family = "serif"),
          plot.title = element_text(size = 35, face = "bold", hjust = 0.5, vjust = 4, margin = margin(30, 0, 0, 0)),
          legend.title = element_blank(),
          legend.background = element_blank(),
          legend.box.background = element_rect(colour = "black"),
          legend.key = element_rect(fill = NA, colour = NA),
          legend.key.width = unit(1.1, "cm"),
          legend.margin = margin(10, 15, 10, 15),
          axis.ticks = element_line(size = 1.5),
          axis.ticks.length = unit(-0.2, "cm"),
          axis.text.x = element_text(margin = margin(10, 0, 0, 0)),  
          axis.text.y = element_text(margin = margin(0, 10, 0, 0)),
          axis.text.y.right = element_text(margin = margin(0, 0, 0, 10))) +
    coord_cartesian(ylim = c(min(df_plot$value)-0.05, 
                             max(df_plot$value)+0.05)) +
    guides(
      color = "none"
    )
} 


### Plot a 3 x 3 Matrix of plots:

main_directory <- dirname(rstudioapi::getSourceEditorContext()$path)

dir.create(gsub(" ", "", paste(main_directory, "/Images")), showWarnings = FALSE)

directory <- gsub(" ", "", paste(main_directory, "/Images/", as.character(Sys.Date())))

dir.create(directory, showWarnings = FALSE)

directory <- gsub(" ", "", paste(directory, "/"))

# Blank spaces

blank_space <- ggplot() + theme_void()

modified_plotlist <- list(blank_space, blank_space, blank_space)
group_size <- 3

for (i in seq(1, length(plotlist), by = group_size)) {
  group <- plotlist[i:min(i + group_size - 1, length(plotlist))]
  modified_plotlist <- c(modified_plotlist, list(blank_space, blank_space, blank_space), group)
}

# Add linear, recessive and dominant:

for (i in 1:3) {
  if (i == 1) title_left = "Recessive"
  else if (i == 2) title_left = "Linear"
  else if (i == 3) title_left = "Dominant"
  
  modified_plotlist[[i]] <-  ggplot() + 
    theme_void() + 
    annotate("text", 
             x = 1, 
             y = 1, 
             label = title_left, 
             angle = 90,             
             size = 16,               
             family = "serif",        
             fontface = "bold"        
    ) + 
    xlim(0.5, 1.5) + 
    ylim(0.5, 1.5)
}

# Reorder

reorder_plotlist_columnwise <- function(modified_plotlist, nrow = 3, ncol = 7) {
  matrix_plots <- matrix(modified_plotlist, nrow = nrow, ncol = ncol)
  as.vector(t(matrix_plots))
}

ordered_plotlist <- reorder_plotlist_columnwise(modified_plotlist, nrow = 3, ncol = 7)


# Plot

plot_grid(plotlist = ordered_plotlist,
          nrow = 3,
          ncol = 7,
          rel_widths = c(0.07, 0.08, 1, 0.07, 1, 0.07, 1)
          )

ggsave(file = paste(paste(directory, "Supplementary_Fig_birth_death_parameters", ".png", sep = ""), sep = ""),  units="px", width=1800*2.1, height=960*2.1, dpi = 60*2.1, bg = "white", limitsize = FALSE)



