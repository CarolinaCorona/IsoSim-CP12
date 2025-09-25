# Check which figures are of importance for the study, as heatmaps were found not to be the optimal multimetric figure, rather the boxplots with CI.

  ##################
 #  Input Fitting #
##################

library(ggplot2)
library(gridExtra)
library(grid)

# Check the proper file path
load("./Run-all-10-irreversible-upper100/Adaptation_all-irreversible-upper100.RData")
#load("./Run-all-10-concentrations-5-90-upper100/Adaptation_all_10_concentration-5-90-upper100.RData")
#load("./Run-all-10-no-concentrations-upper100/Adaptation_all_10_no_concentration-upper100.RData")
#load("./Run-all-10-concentrations-upper100/Adaptation_all_10_concentration-upper100.RData")

# Function to adjust color transparency
transparent_color <- function(color, alpha = 0.5) {
  rgb_col <- col2rgb(color)
  rgb(rgb_col[1]/255, rgb_col[2]/255, rgb_col[3]/255, alpha = alpha)
}  

# Define a sequence of time points from 0 to 90 with step 1
time_points_seq <- seq(0, 90, by = 1)
sim_data2 <- list()
for (points in 1:10){
  sim_data2[[points]] <- time_points_seq
}

# Define a function for plotting experimental data and fitted simulations using ggplot2
plot_fitting_gg <- function(met, exp_data, sd_data, sim_data, col, title_suffix, 
                            use_results_fitting = FALSE, results_fitting = NULL, condition = "WT") {
  
  exp_data_offset <- exp_data + epsilon
  sd_data_offset <- sd_data[[met]] + epsilon
  
  # Create a data frame for experimental data
  exp_df <- data.frame(
    Time = times,
    Value = exp_data_offset,
    SD_low = exp_data_offset - sd_data_offset,
    SD_high = exp_data_offset + sd_data_offset
  )
  
  # Base ggplot for experimental data
  p <- ggplot(exp_df, aes(x = Time, y = Value)) +
    geom_point(color = col, fill = col, shape = 21, size = 2) +
    geom_errorbar(aes(ymin = SD_low, ymax = SD_high), width = 0.2, color = col) +
    ggtitle(paste0(met, title_suffix)) +
    scale_x_continuous(breaks = c(0, 5, 10, 15, 30, 60, 90), limits = c(0, 95)) +  # Specify time points and limit
    #xlab("Time (min)") + ylab("13C-Enrichment") +
    theme_minimal() +
    ylim(-0.05, 0.45) +
    theme(plot.title = element_text(hjust = 0.5),
          axis.title.y = element_blank(),
          axis.text.y = element_text(size = 14),
          axis.text.x = element_text(size = 10),
          axis.title.x = element_blank())
  
  
  
  # Add fitted lines if required
  if (use_results_fitting) {
    for (a in 1:10) {
      fe <- function(t) eval(parse(text = results_fitting[[a]][[condition]][["PG3_1-2-3-M3"]]))  # CHANGE HERE FOR WHEN INPUT IS 2PGA
      fit_df <- data.frame(Time = times, Value = fe(times) + epsilon)
      p <- p + geom_line(data = fit_df, aes(x = Time, y = Value), color = col, linewidth = 1.1)
    }
  } else {
    for (j in 1:10) {
      sim_df <- data.frame(Time = times, Value = sim_data[[j]] + epsilon)
      p <- p + geom_line(data = sim_df, aes(x = Time, y = Value), color = col, linewidth = 1.1)
    }
  }
  
  return(p)
}


# Loop for different metabolites and conditions
plot_list <- list()
for (met in c("2PG", "PEP", "3PG")) {
  if (met == "2PG") {
    plot_list[[paste0(met, "-WT")]] <- plot_fitting_gg(met, data_exp_wt[,2], sd_wt, lapply(plotting_vars, function(x) x[["resF_wt"]][,"PG2_1-2-3"]), transparent_color(col_f[1], 0.45), "-WT")
    plot_list[[paste0(met, "-Δcp12")]] <- plot_fitting_gg(met, data_exp_dcp12[,2], sd_dcp12, lapply(plotting_vars, function(x) x[["resF_dcp12"]][,"PG2_1-2-3"]), transparent_color(col_f[2], 0.45), "-Δcp12")
    plot_list[[paste0(met, "-Δcp12::cp12")]] <- plot_fitting_gg(met, data_exp_dcp12_cp12[,2], sd_dcp12_cp12, lapply(plotting_vars, function(x) x[["resF_dcp12_cp12"]][,"PG2_1-2-3"]), transparent_color(col_f[3], 0.45), "-Δcp12::cp12")
    plot_list[[paste0(met, "-Δcp12::cp12::dCysC")]] <- plot_fitting_gg(met, data_exp_dcp12_cp12dCysC[,2], sd_dcp12_cp12dCysC, lapply(plotting_vars, function(x) x[["resF_dcp12_cp12dCysC"]][,"PG2_1-2-3"]), transparent_color(col_f[4], 0.45), "-Δcp12::cp12::dCysC")
    plot_list[[paste0(met, "-Δcp12::cp12::dCysN")]] <- plot_fitting_gg(met, data_exp_dcp12_cp12dCysN[,2], sd_dcp12_cp12dCysN, lapply(plotting_vars, function(x) x[["resF_dcp12_cp12dCysN"]][,"PG2_1-2-3"]), transparent_color(col_f[5], 0.45), "-Δcp12::cp12::dCysN")
    plot_list[[paste0(met, "-Δcp12::cp12::dCysNC")]] <- plot_fitting_gg(met, data_exp_dcp12_cp12dCysNC[,2], sd_dcp12_cp12dCysNC, lapply(plotting_vars, function(x) x[["resF_dcp12_cp12dCysNC"]][,"PG2_1-2-3"]), transparent_color(col_f[6], 0.45), "-Δcp12::cp12::dCysNC")
  } else if (met == "PEP") {
    plot_list[[paste0(met, "-WT")]] <- plot_fitting_gg(met, data_exp_wt[,3], sd_wt, lapply(plotting_vars, function(x) x[["resF_wt"]][,"PEP_1-2-3"]), transparent_color(col_f[1], 0.45), "-WT")
    plot_list[[paste0(met, "-Δcp12")]] <- plot_fitting_gg(met, data_exp_dcp12[,3], sd_dcp12, lapply(plotting_vars, function(x) x[["resF_dcp12"]][,"PEP_1-2-3"]), transparent_color(col_f[2], 0.45), "-Δcp12")
    plot_list[[paste0(met, "-Δcp12::cp12")]] <- plot_fitting_gg(met, data_exp_dcp12_cp12[,3], sd_dcp12_cp12, lapply(plotting_vars, function(x) x[["resF_dcp12_cp12"]][,"PEP_1-2-3"]), transparent_color(col_f[3], 0.45), "-Δcp12::cp12")
    plot_list[[paste0(met, "-Δcp12::cp12::dCysC")]] <- plot_fitting_gg(met, data_exp_dcp12_cp12dCysC[,3], sd_dcp12_cp12dCysC, lapply(plotting_vars, function(x) x[["resF_dcp12_cp12dCysC"]][,"PEP_1-2-3"]), transparent_color(col_f[4], 0.45), "-Δcp12::cp12::dCysC")
    plot_list[[paste0(met, "-Δcp12::cp12::dCysN")]] <- plot_fitting_gg(met, data_exp_dcp12_cp12dCysN[,3], sd_dcp12_cp12dCysN, lapply(plotting_vars, function(x) x[["resF_dcp12_cp12dCysN"]][,"PEP_1-2-3"]), transparent_color(col_f[5], 0.45), "-Δcp12::cp12::dCysN")
    plot_list[[paste0(met, "-Δcp12::cp12::dCysNC")]] <- plot_fitting_gg(met, data_exp_dcp12_cp12dCysNC[,3], sd_dcp12_cp12dCysNC, lapply(plotting_vars, function(x) x[["resF_dcp12_cp12dCysNC"]][,"PEP_1-2-3"]), transparent_color(col_f[6], 0.45), "-Δcp12::cp12::dCysNC")
  } else if (met == "3PG") {
    plot_list[[paste0(met, " WT")]] <- plot_fitting_gg(met, data_input[,"WTK_PG3_1-2-3"], sd_wt, NULL, transparent_color(col_f[1], 0.45), "-WT", use_results_fitting = TRUE, results_fitting = results_fitting, condition = "WT")
    plot_list[[paste0(met, " Δcp12")]] <- plot_fitting_gg(met, data_input[,"dcp12_PG3_1-2-3"], sd_dcp12, NULL, transparent_color(col_f[2], 0.45), "-Δcp12", use_results_fitting = TRUE, results_fitting = results_fitting, condition = "dcp12")
    plot_list[[paste0(met, " Δcp12::cp12")]] <- plot_fitting_gg(met, data_input[,"dcp12cp12_PG3_1-2-3"], sd_dcp12_cp12, NULL, transparent_color(col_f[3], 0.45), "-Δcp12::cp12", use_results_fitting = TRUE, results_fitting = results_fitting, condition = "dcp12_cp12")
    plot_list[[paste0(met, " Δcp12::cp12::dCysC")]] <- plot_fitting_gg(met, data_input[,"dcp12cp12dCysC_PG3_1-2-3"], sd_dcp12_cp12dCysC, NULL, transparent_color(col_f[4], 0.45), "-Δcp12::cp12::dCysC", use_results_fitting = TRUE, results_fitting = results_fitting, condition = "dcp12_cp12dCysC")
    plot_list[[paste0(met, " Δcp12::cp12::dCysN")]] <- plot_fitting_gg(met, data_input[,"dcp12cp12dCysN_PG3_1-2-3"], sd_dcp12_cp12dCysN, NULL, transparent_color(col_f[5], 0.45), "-Δcp12::cp12::dCysN", use_results_fitting = TRUE, results_fitting = results_fitting, condition = "dcp12_cp12dCysN")
    plot_list[[paste0(met, " Δcp12::cp12::dCysNC")]] <- plot_fitting_gg(met, data_input[,"dcp12cp12dCysNC_PG3_1-2-3"], sd_dcp12_cp12dCysNC, NULL, transparent_color(col_f[6], 0.45), "-Δcp12::cp12::dCysNC", use_results_fitting = TRUE, results_fitting = results_fitting, condition = "dcp12_cp12dCysNC")
  }
}


plot_list_3PG <- plot_list[grep("3PG", names(plot_list))]


grid.arrange(
  grobs = plot_list_3PG,
  left = textGrob("13C-Enrichment", rot = 90, gp = gpar(fontsize = 14)),
  bottom = textGrob("Time (min)", gp = gpar(fontsize = 14)),
  ncol = length(plot_list_3PG)  # Horizontal arrangement (or use nrow for vertical)
)




# Plot for WT, Δcp12, and additional species for 3PG and 2PG
plot_list <- list()
for (met in c("2PG", "3PG")) {
  if (met == "2PG") {
    # Plot WT for PEP using plotting_vars with transparency
    plot_list[[paste0(met, "-WT")]] <- plot_fitting_gg(met, data_exp_wt[,2], sd_wt, lapply(plotting_vars, function(x) x[["resF_wt"]]), transparent_color(col_f[1], 0.45), "-WT")
    
    # Plot Δcp12 for PEP using plotting_vars with transparency
    plot_list[[paste0(met, "-Δcp12")]] <- plot_fitting_gg(met, data_exp_dcp12[,2], sd_dcp12, lapply(plotting_vars, function(x) x[["resF_dcp12"]]), transparent_color(col_f[2], 0.45), "-Δcp12")
    
    # Plot Δcp12::cp12 for 2PG using plotting_vars with transparency
    plot_list[[paste0(met, "-Δcp12::cp12")]] <- plot_fitting_gg(met, data_exp_dcp12_cp12[,2], sd_dcp12_cp12, lapply(plotting_vars, function(x) x[["resF_dcp12_cp12"]]), transparent_color(col_f[3], 0.45), "-Δcp12::cp12")
    
    # Additional species with transparency
    plot_list[[paste0(met, "-Δcp12::cp12::dCysC")]] <- plot_fitting_gg(met, data_exp_dcp12_cp12dCysC[,2], sd_dcp12_cp12dCysC, lapply(plotting_vars, function(x) x[["resF_dcp12_cp12dCysC"]]), transparent_color(col_f[4], 0.45), "-Δcp12::cp12::dCysC")
    plot_list[[paste0(met, "-Δcp12::cp12::dCysN")]] <-  plot_fitting_gg(met, data_exp_dcp12_cp12dCysN[,2], sd_dcp12_cp12dCysN, lapply(plotting_vars, function(x) x[["resF_dcp12_cp12dCysN"]]), transparent_color(col_f[5], 0.45), "-Δcp12::cp12::dCysN")
    plot_list[[paste0(met, "-Δcp12::cp12::dCysNC")]] <- plot_fitting_gg(met, data_exp_dcp12_cp12dCysNC[,2], sd_dcp12_cp12dCysNC, lapply(plotting_vars, function(x) x[["resF_dcp12_cp12dCysNC"]]), transparent_color(col_f[6], 0.45), "-Δcp12::cp12::dCysNC")
    
  } else if (met == "3PG") {
    # Plot WT for 2PG using results_fitting with transparency
    plot_list[[paste0(met, "-WT")]] <- plot_fitting_gg(met, data_input[,"WTK_PG3_1-2-3"], sd_wt, results_fitting, transparent_color(col_f[1], 0.45), "-WT", use_results_fitting = TRUE, results_fitting = results_fitting, condition = "WT")
    
    # Plot Δcp12 for 2PG using results_fitting with transparency
    plot_list[[paste0(met, "-Δcp12")]] <- plot_fitting_gg(met, data_input[,"dcp12_PG3_1-2-3"], sd_dcp12, results_fitting, transparent_color(col_f[2], 0.45), "-Δcp12", use_results_fitting = TRUE, results_fitting = results_fitting, condition = "dcp12")
    
    # Plot Δcp12::cp12 for 2PG using results_fitting with transparency
    plot_list[[paste0(met, "-Δcp12::cp12")]] <- plot_fitting_gg(met, data_input[,"dcp12cp12_PG3_1-2-3"], sd_dcp12_cp12, results_fitting, transparent_color(col_f[3], 0.45), "-Δcp12::cp12", use_results_fitting = TRUE, results_fitting = results_fitting, condition = "dcp12_cp12")
    
    # Additional species with transparency
    plot_list[[paste0(met, "-Δcp12::cp12::dCysC")]] <- plot_fitting_gg(met, data_input[,"dcp12cp12dCysC_PG3_1-2-3"], sd_dcp12_cp12dCysC, results_fitting, transparent_color(col_f[4], 0.45), "-Δcp12::cp12::dCysC", use_results_fitting = TRUE, results_fitting = results_fitting, condition = "dcp12_cp12dCysC")
    plot_list[[paste0(met, "-Δcp12::cp12::dCysN")]] <- plot_fitting_gg(met, data_input[,"dcp12cp12dCysN_PG3_1-2-3"], sd_dcp12_cp12dCysN, results_fitting, transparent_color(col_f[5], 0.45), "-Δcp12::cp12::dCysN", use_results_fitting = TRUE, results_fitting = results_fitting, condition = "dcp12_cp12dCysN")
    plot_list[[paste0(met, "-Δcp12::cp12::dCysNC")]] <- plot_fitting_gg(met, data_input[,"dcp12cp12dCysNC_PG3_1-2-3"], sd_dcp12_cp12dCysNC, results_fitting, transparent_color(col_f[6], 0.45), "-Δcp12::cp12::dCysNC", use_results_fitting = TRUE, results_fitting = results_fitting, condition = "dcp12_cp12dCysNC")
  }
}
plot_list_3PG <- plot_list[grep("3PG", names(plot_list))]


grid.arrange(
  grobs = plot_list_3PG,
  left = textGrob("13C-Enrichment", rot = 90, gp = gpar(fontsize = 14)),
  bottom = textGrob("Time (min)", gp = gpar(fontsize = 14)),
  ncol = length(plot_list_3PG)  # Horizontal arrangement (or use nrow for vertical)
)




# Plot for WT, Δcp12, and additional species for PEP and 2PG
plot_list <- list()
for (met in c("2PG", "PEP")) {
  if (met == "PEP") {
    # Plot WT for PEP using plotting_vars with transparency
    plot_list[[paste0(met, "-WT")]] <- plot_fitting_gg(met, data_exp_wt[,2], sd_wt, lapply(plotting_vars, function(x) x[["resF_wt"]]), transparent_color(col_f[1], 0.45), "-WT")

    # Plot Δcp12 for PEP using plotting_vars with transparency
    plot_list[[paste0(met, "-Δcp12")]] <- plot_fitting_gg(met, data_exp_dcp12[,2], sd_dcp12, lapply(plotting_vars, function(x) x[["resF_dcp12"]]), transparent_color(col_f[2], 0.45), "-Δcp12")

    # Plot Δcp12::cp12 for 2PG using plotting_vars with transparency
    plot_list[[paste0(met, "-Δcp12::cp12")]] <- plot_fitting_gg(met, data_exp_dcp12_cp12[,2], sd_dcp12_cp12, lapply(plotting_vars, function(x) x[["resF_dcp12_cp12"]]), transparent_color(col_f[3], 0.45), "-Δcp12::cp12")

    # Additional species with transparency
    plot_list[[paste0(met, "-Δcp12::cp12::dCysC")]] <- plot_fitting_gg(met, data_exp_dcp12_cp12dCysC[,2], sd_dcp12_cp12dCysC, lapply(plotting_vars, function(x) x[["resF_dcp12_cp12dCysC"]]), transparent_color(col_f[4], 0.45), "-Δcp12::cp12::dCysC")
    plot_list[[paste0(met, "-Δcp12::cp12::dCysN")]] <-  plot_fitting_gg(met, data_exp_dcp12_cp12dCysN[,2], sd_dcp12_cp12dCysN, lapply(plotting_vars, function(x) x[["resF_dcp12_cp12dCysN"]]), transparent_color(col_f[5], 0.45), "-Δcp12::cp12::dCysN")
    plot_list[[paste0(met, "-Δcp12::cp12::dCysNC")]] <- plot_fitting_gg(met, data_exp_dcp12_cp12dCysNC[,2], sd_dcp12_cp12dCysNC, lapply(plotting_vars, function(x) x[["resF_dcp12_cp12dCysNC"]]), transparent_color(col_f[6], 0.45), "-Δcp12::cp12::dCysNC")

  } else if (met == "2PG") {
    # Plot WT for 2PG using results_fitting with transparency
    plot_list[[paste0(met, "-WT")]] <- plot_fitting_gg(met, data_input[,"WTK_PG2_1-2-3"], sd_wt, results_fitting, transparent_color(col_f[1], 0.45), "-WT", use_results_fitting = TRUE, results_fitting = results_fitting, condition = "WT")
    
    # Plot Δcp12 for 2PG using results_fitting with transparency
    plot_list[[paste0(met, "-Δcp12")]] <- plot_fitting_gg(met, data_input[,"dcp12_PG2_1-2-3"], sd_dcp12, results_fitting, transparent_color(col_f[2], 0.45), "-Δcp12", use_results_fitting = TRUE, results_fitting = results_fitting, condition = "dcp12")
    
    # Plot Δcp12::cp12 for 2PG using results_fitting with transparency
    plot_list[[paste0(met, "-Δcp12::cp12")]] <- plot_fitting_gg(met, data_input[,"dcp12cp12_PG2_1-2-3"], sd_dcp12_cp12, results_fitting, transparent_color(col_f[3], 0.45), "-Δcp12::cp12", use_results_fitting = TRUE, results_fitting = results_fitting, condition = "dcp12_cp12")
    
    # Additional species with transparency
    plot_list[[paste0(met, "-Δcp12::cp12::dCysC")]] <- plot_fitting_gg(met, data_input[,"dcp12cp12dCysC_PG2_1-2-3"], sd_dcp12_cp12dCysC, results_fitting, transparent_color(col_f[4], 0.45), "-Δcp12::cp12::dCysC", use_results_fitting = TRUE, results_fitting = results_fitting, condition = "dcp12_cp12dCysC")
    plot_list[[paste0(met, "-Δcp12::cp12::dCysN")]] <- plot_fitting_gg(met, data_input[,"dcp12cp12dCysN_PG2_1-2-3"], sd_dcp12_cp12dCysN, results_fitting, transparent_color(col_f[5], 0.45), "-Δcp12::cp12::dCysN", use_results_fitting = TRUE, results_fitting = results_fitting, condition = "dcp12_cp12dCysN")
    plot_list[[paste0(met, "-Δcp12::cp12::dCysNC")]] <- plot_fitting_gg(met, data_input[,"dcp12cp12dCysNC_PG2_1-2-3"], sd_dcp12_cp12dCysNC, results_fitting, transparent_color(col_f[6], 0.45), "-Δcp12::cp12::dCysNC", use_results_fitting = TRUE, results_fitting = results_fitting, condition = "dcp12_cp12dCysNC")
  }
}

plot_list_2PG <- plot_list[grep("2PG", names(plot_list))]


grid.arrange(
  grobs = plot_list_2PG,
  left = textGrob("13C-Enrichment", rot = 90, gp = gpar(fontsize = 14)),
  bottom = textGrob("Time (min)", gp = gpar(fontsize = 14)),
  ncol = length(plot_list_2PG)  # Horizontal arrangement (or use nrow for vertical)
)



  ###################
 # Goodness of fit #
##################

# Load necessary libraries
library(ggplot2)
library(reshape2)
library(dplyr)

# Assuming results_flux and result_vars2 are pre-defined and accessible

# Reshape the flux data into a long format
df_flux <- melt(as.data.frame(results_flux), variable.name = "Strain", value.name = "Flux")

# Extract the chi-square conclusion for each run of each strain
fit_conclusions_raw <- c(
  sapply(1:10, function(i) result_vars2[[i]][[1]]$chi2$conclusion),  # WT
  sapply(1:10, function(i) result_vars2[[i]][[2]]$chi2$conclusion),  # dcp12
  sapply(1:10, function(i) result_vars2[[i]][[3]]$chi2$conclusion),  # dcp12::cp12
  sapply(1:10, function(i) result_vars2[[i]][[4]]$chi2$conclusion),  # dcp12::cp12dCysC
  sapply(1:10, function(i) result_vars2[[i]][[5]]$chi2$conclusion),  # dcp12::cp12dCysN
  sapply(1:10, function(i) result_vars2[[i]][[6]]$chi2$conclusion)   # dcp12::cp12dCysNC
)

# Convert detailed fit conclusions to "Yes" or "No"
fit_conclusions <- ifelse(grepl("fits the data good enough", fit_conclusions_raw), "Yes", "No")

# Results of chi-square for each strain and run
chi2_value <- c(
  sapply(1:10, function(i) result_vars2[[i]][[1]]$chi2$"khi2 value"),  # WT
  sapply(1:10, function(i) result_vars2[[i]][[2]]$chi2$"khi2 value"),  # dcp12
  sapply(1:10, function(i) result_vars2[[i]][[3]]$chi2$"khi2 value"),  # dcp12::cp12
  sapply(1:10, function(i) result_vars2[[i]][[4]]$chi2$"khi2 value"),  # dcp12::cp12dCysC
  sapply(1:10, function(i) result_vars2[[i]][[5]]$chi2$"khi2 value"),  # dcp12::cp12dCysN
  sapply(1:10, function(i) result_vars2[[i]][[6]]$chi2$"khi2 value")   # dcp12::cp12dCysNC
)

chi2_reduced_value <- c(
  sapply(1:10, function(i) result_vars2[[i]][[1]]$chi2$"khi2 reduced value"),  # WT
  sapply(1:10, function(i) result_vars2[[i]][[2]]$chi2$"khi2 reduced value"),  # dcp12
  sapply(1:10, function(i) result_vars2[[i]][[3]]$chi2$"khi2 reduced value"),  # dcp12::cp12
  sapply(1:10, function(i) result_vars2[[i]][[4]]$chi2$"khi2 reduced value"),  # dcp12::cp12dCysC
  sapply(1:10, function(i) result_vars2[[i]][[5]]$chi2$"khi2 reduced value"),  # dcp12::cp12dCysN
  sapply(1:10, function(i) result_vars2[[i]][[6]]$chi2$"khi2 reduced value")   # dcp12::cp12dCysNC
)

p_value <- c(
  sapply(1:10, function(i) result_vars2[[i]][[1]]$chi2$"p-value, i.e. P(X^2<=value)"),  # WT
  sapply(1:10, function(i) result_vars2[[i]][[2]]$chi2$"p-value, i.e. P(X^2<=value)"),  # dcp12
  sapply(1:10, function(i) result_vars2[[i]][[3]]$chi2$"p-value, i.e. P(X^2<=value)"),  # dcp12::cp12
  sapply(1:10, function(i) result_vars2[[i]][[4]]$chi2$"p-value, i.e. P(X^2<=value)"),  # dcp12::cp12dCysC
  sapply(1:10, function(i) result_vars2[[i]][[5]]$chi2$"p-value, i.e. P(X^2<=value)"),  # dcp12::cp12dCysN
  sapply(1:10, function(i) result_vars2[[i]][[6]]$chi2$"p-value, i.e. P(X^2<=value)")   # dcp12::cp12dCysNC
)

degreesf <- c(
  sapply(1:10, function(i) result_vars2[[i]][[1]]$chi2$"degrees of freedom"),  # WT
  sapply(1:10, function(i) result_vars2[[i]][[2]]$chi2$"degrees of freedom"),  # dcp12
  sapply(1:10, function(i) result_vars2[[i]][[3]]$chi2$"degrees of freedom"),  # dcp12::cp12
  sapply(1:10, function(i) result_vars2[[i]][[4]]$chi2$"degrees of freedom"),  # dcp12::cp12dCysC
  sapply(1:10, function(i) result_vars2[[i]][[5]]$chi2$"degrees of freedom"),  # dcp12::cp12dCysN
  sapply(1:10, function(i) result_vars2[[i]][[6]]$chi2$"degrees of freedom")   # dcp12::cp12dCysNC
)

# Create a simplified data frame of chi-square results for each condition
chi2_data <- data.frame(
  Strain = rep(c("WT", "dcp12", "dcp12::cp12", "dcp12::cp12dCysC", "dcp12::cp12dCysN", "dcp12::cp12dCysNC"), each = 10),
  Run = rep(1:10, times = 6),
  Chi2 = chi2_value,  
  Chi2_Reduced = chi2_reduced_value,  
  P_Value = p_value 
)

# Merge chi2_data with the main flux data
df_flux <- cbind(df_flux, chi2_data[, c("Run","Chi2", "Chi2_Reduced", "P_Value")], fit_conclusions)
df_flux[,5:7]
df_flux$Chi2_Reduced

max(df_flux$Chi2_Reduced)

order(df_flux$Chi2_Reduced)

ggplot(df_flux, aes(x = Run, y = Chi2_Reduced, color = fit_conclusions)) +
  geom_point(size = 3) +
  geom_line(aes(group = Strain), alpha = 0.5) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "#0072B2") +
  scale_x_continuous(breaks = c(1:10)) +
  ylim(-0.05, 3) +
  facet_wrap(~Strain, scales = "free_x") +
  labs(
    title = "Reduced Chi-Square (χ²/df) by Strain and Run",
    x = "Run",
    y = "Reduced Chi-Square (χ²/df)",
    color = "Is the fit good?"
  ) +
  scale_color_manual(values = c("Yes" = "#009E73", "No" = "#D55E00")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10), 
        plot.title = element_blank(), axis.text.y = element_text(size = 12),
        legend.position = "bottom")


# library(writexl)
# write_xlsx(df_flux, "./Run-2pg-pep-10-irreversible-upper100/2pg-pep-irreversible-upper100.xlsx")
# write_xlsx(as.data.frame(results_flux), "./Run-2pg-pep-10-irreversible-upper100/Result fluxes 2pg pep irreversible new upper 100.xlsx")


  ###################
 # Estimated fluxes#
###################

col_factor <- c("WT" = "#1f77b4", 
                "dcp12" = "#ff7f0e", 
                "dcp12::cp12" = "#2ca02c", 
                "dcp12::cp12dCysC" = "#d62728", 
                "dcp12::cp12dCysN" = "#9467bd", 
                "dcp12::cp12dCysNC" = "#8c564b")

# Define a small tolerance for floating-point comparison
tolerance <- 1e-10

# Define a function to calculate the Flux Log
calculate_flux_log <- function(flux) {
  if (abs(flux - 0.0001) < tolerance) {
    return(log10(2*flux))
  } else {
    return(log10(flux))
  }
}

# Apply the function to the Flux column to create the new Flux Log column
df_flux <- df_flux %>%
  mutate(Flux_Log = sapply(Flux, calculate_flux_log))

# Display the updated data frame
print(df_flux)

ggplot(df_flux, aes(x = Strain, y = Flux_Log, fill = Strain)) +
  geom_boxplot(alpha = 0.75, outlier.shape = 16, outlier.size = 2) +
  scale_fill_manual(values = col_factor) +
  theme_minimal() +
  labs(title = "",
       x = "Strain",
       y = "LOG10(Flux)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5))+
  ylim(-4, 2.1)


# Create the plot with position_jitter() and white outlines
flux_log <- ggplot(df_flux, aes(x = Strain, y = Flux_Log)) +
  geom_point(
    aes(fill = Strain),  # Use fill for the interior color
    position = position_jitter(width = 0.4, height = 0),  # Control jittering
    alpha = 0.7,  # Adjust transparency
    shape = 21,   # Shape with both fill and outline (21-25)
    size = 3,     # Size of the points
    color = "white"  # Outline color (white)
  ) +
  scale_fill_manual(values = col_factor) +  # Use your custom color palette for fill
  scale_y_continuous(breaks = seq(-4, 2, by = 1), limits = c(-4, 2.1)) +
  labs(
    title = "", 
    x = "Strain", 
    y = "LOG10(Flux)"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8), 
    plot.title = element_text(hjust = 0.5),
    legend.position = "none"
  )

# Print the plot
print(flux_log)


  #######################
 # Confidence Interval #
######################

# Load necessary libraries
library(ggplot2)
library(reshape2)
library(dplyr)

load("./Run-all-10-no-concentrations-upper100/Adaptation_all_10_no_concentration-upper100.RData")
#load("./Run-3pg-2pg-10-no-concentrations-upper100/Adaptation-3pg-2pg-10-no-concentrations-upper100.RData")
#load("./Run-2pg-pep-10-no-concentrations-upper100/Adaptation-2pg-pep-10-no-concentrations-upper100.RData")

#load("./Run-all-10-concentrations-upper100/Adaptation_all_10_concentration-upper100.RData")
#load("./Run-3pg-2pg-10-concentrations-upper100/Adaptation_3pg_2pg_10_concentration-upper100.RData")

#load("./Run-all-10-concentrations-5-90-upper100/Adaptation_all_10_concentration-5-90-upper100.RData")
#load("./Run-3pg-2pg-10-concentrations-5-90-upper100/Adaptation_3pg_2pg_10_concentration-5-90-upper100.RData")

#load("./Run-all-10-irreversible-upper100/Adaptation_all-irreversible-upper100.RData")
#load("./Run-3pg-2pg-10-irreversible-upper100/Adaptation_3pg-2pg-irreversible-upper100.RData")
#load("./Run-2pg-pep-10-irreversible-upper100/Adaptation_2pg-pep-irreversible-upper100.RData")


#Gotta do a 3D matrix to store some info for plotting
# Define the strains and metrics
strains <- c("WT", "dcp12", "dcp12::cp12", "dcp12::cp12dCysC", "dcp12::cp12dCysN", "dcp12::cp12dCysNC")
metrics <- c("opt", "mean", "median", "sd", "ci_2.5", "ci_97.5")

# Initialize a 3D array to store results for all 10 runs
# Dimensions: strains x metrics x runs
conf_array <- array(NA, dim = c(length(strains), length(metrics), 10),
                    dimnames = list(strains, metrics, paste0("Run_", 1:10)))

# Loop through the 10 runs and fill the array
for (r in 1:10) {
  for (i in 1:length(strains)) {
    strain <- strains[i]
    
    # Extract the results for the current strain and run
    if (strain == "WT") {
      conf_array[i, , r] <- result_vars2[[r]]$resF_wt$sens$summary["v1", ]
    } else if (strain == "dcp12") {
      conf_array[i, , r] <- result_vars2[[r]]$resF_dcp12$sens$summary["v1", ]
    } else if (strain == "dcp12::cp12") {
      conf_array[i, , r] <- result_vars2[[r]]$resF_dcp12_cp12$sens$summary["v1", ]
    } else if (strain == "dcp12::cp12dCysC") {
      conf_array[i, , r] <- result_vars2[[r]]$resF_dcp12_cp12dCysC$sens$summary["v1", ]
    } else if (strain == "dcp12::cp12dCysN") {
      conf_array[i, , r] <- result_vars2[[r]]$resF_dcp12_cp12dCysN$sens$summary["v1", ]
    } else if (strain == "dcp12::cp12dCysNC") {
      conf_array[i, , r] <- result_vars2[[r]]$resF_dcp12_cp12dCysNC$sens$summary["v1", ]
    }
  }
}

# # Check the result for a specific run (e.g., Run 1)
# print(conf_array[, , 1])
# 
# # Check the result for a specific strain across all runs (e.g., WT)
# print(conf_array["dcp12", , ])
# 
# # Check the result for a specific metric across all runs (e.g., mean)
# print(conf_array[, "mean", ])



# Extract data for plotting
plotci_data <- data.frame(
  Strain = character(),
  Run = numeric(),
  Median = numeric(),
  Mean = numeric(),
  CI_lower = numeric(),
  CI_upper = numeric()
)

# Loop through the 10 runs and strains to extract data
for (r in 1:10) {
  for (i in 1:length(strains)) {
    strain <- strains[i]
    
    # Extract the mean, ci_2.5, and ci_97.5 for the current strain and run
    if (strain == "WT") {
      mean_val <- result_vars2[[r]]$resF_wt$sens$summary["v1", "mean"]
      ci_lower <- result_vars2[[r]]$resF_wt$sens$summary["v1", "ci_2.5"]
      ci_upper <- result_vars2[[r]]$resF_wt$sens$summary["v1", "ci_97.5"]
    } else if (strain == "dcp12") {
      mean_val <- result_vars2[[r]]$resF_dcp12$sens$summary["v1", "mean"]
      ci_lower <- result_vars2[[r]]$resF_dcp12$sens$summary["v1", "ci_2.5"]
      ci_upper <- result_vars2[[r]]$resF_dcp12$sens$summary["v1", "ci_97.5"]
    } else if (strain == "dcp12::cp12") {
      mean_val <- result_vars2[[r]]$resF_dcp12_cp12$sens$summary["v1", "mean"]
      ci_lower <- result_vars2[[r]]$resF_dcp12_cp12$sens$summary["v1", "ci_2.5"]
      ci_upper <- result_vars2[[r]]$resF_dcp12_cp12$sens$summary["v1", "ci_97.5"]
    } else if (strain == "dcp12::cp12dCysC") {
      mean_val <- result_vars2[[r]]$resF_dcp12_cp12dCysC$sens$summary["v1", "mean"]
      ci_lower <- result_vars2[[r]]$resF_dcp12_cp12dCysC$sens$summary["v1", "ci_2.5"]
      ci_upper <- result_vars2[[r]]$resF_dcp12_cp12dCysC$sens$summary["v1", "ci_97.5"]
    } else if (strain == "dcp12::cp12dCysN") {
      mean_val <- result_vars2[[r]]$resF_dcp12_cp12dCysN$sens$summary["v1", "mean"]
      ci_lower <- result_vars2[[r]]$resF_dcp12_cp12dCysN$sens$summary["v1", "ci_2.5"]
      ci_upper <- result_vars2[[r]]$resF_dcp12_cp12dCysN$sens$summary["v1", "ci_97.5"]
    } else if (strain == "dcp12::cp12dCysNC") {
      mean_val <- result_vars2[[r]]$resF_dcp12_cp12dCysNC$sens$summary["v1", "mean"]
      ci_lower <- result_vars2[[r]]$resF_dcp12_cp12dCysNC$sens$summary["v1", "ci_2.5"]
      ci_upper <- result_vars2[[r]]$resF_dcp12_cp12dCysNC$sens$summary["v1", "ci_97.5"]
    }
    
    # Add the data to the plot_data data frame
    plotci_data <- rbind(plotci_data, data.frame(
      Strain = strain,
      Run = r,
      Mean = mean_val,
      CI_lower = ci_lower,
      CI_upper = ci_upper
    ))
  }
}


# Calculate CI width and add as new column
plotci_data$CI_width <- plotci_data$CI_upper - plotci_data$CI_lower

# Add yes/no column (50% threshold of Mean flux)
plotci_data$CI_acceptable <- ifelse(plotci_data$CI_width < 0.5 * plotci_data$Mean, "yes", "no")

# View the modified dataframe
plotci_data


library(ggplot2)
library(ggtext) # Needed for element_markdown() to render bold text

# Create a new column with formatted CI text (bold if acceptable)
plotci_data <- plotci_data %>%
  mutate(
    CI_formatted = ifelse(
      CI_acceptable == "yes",
      sprintf("<b>%.3f</b>", CI_width), # Bold for acceptable
      sprintf("%.3f", CI_width)         # Regular for others
    )
  )

plotci_data



############################# ENHANCED HEATMAP ######################################

library(ggplot2)
library(ggtext)
library(dplyr)
library(readxl)

merged_flux <- read_excel("./Merged_Max_Min_Flux.xlsx")


# Step 1: Calculate the true lowest min_flux across all metabolites
merged_flux <- merged_flux %>%
  rowwise() %>%  # Ensure calculations are done row-by-row
  mutate(
    lowest_min_flux = min(min_flux_3PGA, min_flux_2PGA, min_flux_PEP, na.rm = TRUE)
  ) %>%
  ungroup()

# Step 2: Merge and flag - fixing NA issue for WT
# Step 2: Merge with plot data and flag
plotci_data <- plotci_data %>%
  left_join(
    merged_flux %>% select(Strain, lowest_min_flux),
    by = "Strain"
  ) %>%
  mutate(
    flag = ifelse(is.na(lowest_min_flux), FALSE, Mean < lowest_min_flux)
  )

# Step 3: Format labels with * for flagged cells
plotci_data <- plotci_data %>%
  mutate(
    Run = factor(Run, levels = 1:10),  # Force Run as ordered factor
    CI_formatted = case_when(
      is.na(Mean) ~ NA_character_,
      CI_acceptable == "yes" & flag ~ sprintf("<b>%.3f*</b>", CI_width),
      CI_acceptable == "yes" ~ sprintf("<b>%.3f</b>", CI_width),
      flag ~ sprintf("%.3f*", CI_width),
      TRUE ~ sprintf("%.3f", CI_width)
    )
  )

# Step 4: Plot with all fixes
NEW_heatmap <- function(data, title) {
  # Calculate normalization factors
  norm_factors <- data %>%
    group_by(Strain) %>%
    summarise(max_val = max(Mean, na.rm = TRUE))
  
  # Apply normalization
  plot_data <- data %>%
    left_join(norm_factors, by = "Strain") %>%
    mutate(Normalized_Mean = Mean / max_val)
  
  p <- ggplot(plot_data, aes(x = Run, y = Strain, fill = Normalized_Mean)) +
    geom_tile(color = "white", linewidth = 0.3) +
    geom_richtext(
      aes(label = CI_formatted),
      size = 3,
      fill = NA, 
      label.color = NA,
      color = "black"
    ) +
    scale_fill_gradientn(
      colours = c("antiquewhite1", "darkorchid3"),
      limits = c(0, 1),
      labels = scales::percent_format(),
      na.value = "gray90"
    ) +
    scale_x_discrete(
      breaks = 1:10,
      expand = c(0, 0)
    ) +
    labs(
      x = "",
      y = "Strain",
      fill = "Relative Mean Flux\n(% of max)",
      title = title
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      axis.text.y = element_text(face = "bold"),
      panel.grid = element_blank(),
      plot.title = element_text(size = 10),
      legend.position = "right"
    ) +
    coord_fixed(ratio = 0.8)
  
  print(p)
}

#p_4_4_4 <- plotci_data
#save.image("./p_4_4_4.RData")


# Load all your saved .RData files (adjust paths as needed)
load("./p_1.RData")
load("./p_1_1.RData")
load("./p_1_1_1.RData")
load("./p_2.RData")
load("./p_2_2.RData")
load("./p_3.RData")
load("./p_3_3.RData")
load("./p_4.RData")
load("./p_4_4.RData")
load("./p_4_4_4.RData")

library(ggplot2)
library(ggtext)
library(patchwork) # For combining plots
library(dplyr)

# Create plots (only left column shows strain names)
p1 <- NEW_heatmap(p_1, "3PGA ↔ 2PGA ↔ PEP → ∅ No Conc 0-90 min.")
p1_1 <- NEW_heatmap(p_1_1, "3PGA ↔ 2PGA → ∅ No Conc 0-90 min.")
p1_1_1 <- NEW_heatmap(p_1_1_1, "2PGA ↔ PEP → ∅ No Conc 0-90 min.")

p2 <- NEW_heatmap(p_2, "3PGA ↔ 2PGA ↔ PEP → ∅ Conc 0-90 min.")
p2_2 <- NEW_heatmap(p_2_2, "3PGA ↔ 2PGA → ∅ Conc 0-90 min.")

p3 <- NEW_heatmap(p_3, "3PGA ↔ 2PGA ↔ PEP → ∅ Conc 5-90 min.")
p3_3 <- NEW_heatmap(p_3_3, "3PGA ↔ 2PGA → ∅ Conc 5-90 min.")

p4 <- NEW_heatmap(p_4, "3PGA → 2PGA → PEP → ∅ Conc 0-90 min.")
p4_4 <- NEW_heatmap(p_4_4, "3PGA → 2PGA → ∅ Conc 0-90 min")
p4_4_4 <- NEW_heatmap(p_4_4_4, "2PGA → PEP → ∅ Conc 0-90 min")


# No Conc 0-90

Scenario1 <- ( p1 / p1_1 /  p1_1_1 ) +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "right",
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 9),
    plot.margin = margin(0, 0, 0, 0, unit = "pt")
  )

# Conc 0-90
Scenario2 <- ( p2 / p2_2 ) +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "right",
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 9),
    plot.margin = margin(0, 0, 0, 0, unit = "pt")
  )  

# Conc 5-90
Scenario3 <- ( p3 / p3_3 ) +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "right",
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 9),
    plot.margin = margin(0, 0, 0, 0, unit = "pt")
  ) 

# Conc irrev 0-90
Scenario4 <- ( p4 / p4_4 /  p4_4_4 ) +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "right",
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 9),
    plot.margin = margin(0, 0, 0, 0, unit = "pt")
  )



save.image("NewHeatmaps.RData")




###################### PLOT TOGETHER #########################

library(ggplot2)
library(ggtext)
library(patchwork) # For combining plots
library(dplyr)

# Load all your saved .RData files (adjust paths as needed)
load("./plotci_1.RData")
load("./plotci_1_1.RData")
load("./plotci_1_1_1.RData")
load("./plotci_2.RData")
load("./plotci_2_2.RData")
load("./plotci_3.RData")
load("./plotci_3_3.RData")
load("./plotci_4.RData")
load("./plotci_4_4.RData")
load("./plotci_4_4_4.RData")



# List all plotci_data objects
all_data <- list(
  plotci_data_1, plotci_data_1_1, plotci_data_1_1_1,
  plotci_data_2, plotci_data_2_2,
  plotci_data_3, plotci_data_3_3,
  plotci_data_4, plotci_data_4_4, plotci_data_4_4_4
)


# Find global min and max of Mean values
global_min <- min(sapply(all_data, function(x) min(x$Mean, na.rm = TRUE)))
global_max <- max(sapply(all_data, function(x) max(x$Mean, na.rm = TRUE)))



library(patchwork)

# Modified heatmap function with axis control
create_heatmap <- function(data, title, show_strain_names = TRUE) {
  p <- ggplot(data, aes(x = Run, y = Strain, fill = Mean)) +
    geom_tile() +
    geom_richtext(
      aes(label = CI_formatted),
      size = 2.5,
      fill = NA, label.color = NA,
      color = "black"
    ) +
    scale_fill_gradient(
      low = "orange", high = "purple2",
      limits = c(global_min, global_max)
    ) +
    scale_x_continuous(breaks = 1:10, limits = c(0, 11)) +
    labs(title = title, x = element_blank(), y = if(show_strain_names) "Strain" else NULL) +
    #theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
      axis.text.y = if(show_strain_names) {
        element_text(size = 7, face = "bold")
      } else {
        element_blank()
      },
      plot.title = element_text(hjust = 0.5, size = 8),
      legend.position = "none",
      plot.margin = margin(2, 2, -5, 2)
    )
  return(p)
}

# Create plots (only left column shows strain names)
p1 <- create_heatmap(plotci_data_1, "3PGA ↔ 2PGA ↔ PEP → ∅ No Conc 0-90 min.", show_strain_names = TRUE)
p1_1 <- create_heatmap(plotci_data_1_1, "3PGA ↔ 2PGA → ∅ No Conc 0-90 min.", T)
p1_1_1 <- create_heatmap(plotci_data_1_1_1, "2PGA ↔ PEP → ∅ No Conc 0-90 min.", T)

p2 <- create_heatmap(plotci_data_2, "3PGA ↔ 2PGA ↔ PEP → ∅ Conc 0-90 min.", TRUE)
p2_2 <- create_heatmap(plotci_data_2_2, "3PGA ↔ 2PGA → ∅ Conc 0-90 min.", T)

p3 <- create_heatmap(plotci_data_3, "3PGA ↔ 2PGA ↔ PEP → ∅ Conc 5-90 min.", TRUE)
p3_3 <- create_heatmap(plotci_data_3_3, "3PGA ↔ 2PGA → ∅ Conc 5-90 min.", T)

p4 <- create_heatmap(plotci_data_4, "3PGA → 2PGA → PEP → ∅ Conc 0-90 min.", TRUE)
p4_4 <- create_heatmap(plotci_data_4_4, "3PGA → 2PGA → ∅ Conc 0-90 min", T)
p4_4_4 <- create_heatmap(plotci_data_4_4_4, "2PGA → PEP → ∅ Conc 0-90 min", T)


# No Conc 0-90

Scenario1 <- ( p1 / p1_1 /  p1_1_1 ) +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "right",
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 9)
  )

# Conc 0-90
Scenario2 <- ( p2 / p2_2 ) +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "right",
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 9)
  )  

# Conc 5-90
Scenario3 <- ( p3 / p3_3 ) +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "right",
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 9)
  ) 

# Conc irrev 0-90
Scenario4 <- ( p4 / p4_4 /  p4_4_4 ) +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "right",
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 9)
  )

# Combine with negative spacing between plots
combined_plot <- (
  (p1 + theme(plot.margin = margin(0, -25, -15, 0))) |
    (p1_1 + theme(plot.margin = margin(0, -25, -15, 0))) |
    p1_1_1
) / (
  (p2 + theme(plot.margin = margin(-3, -25, -15, 0))) |
    (p2_2 + theme(plot.margin = margin(-3, -25, -15, 0))) |
    plot_spacer()
) / (
  (p3 + theme(plot.margin = margin(-3, -25, -15, 0))) |
    (p3_3 + theme(plot.margin = margin(-3, -25, -15, 0))) |
    plot_spacer()
) / (
  (p4 + theme(plot.margin = margin(-3, -25, -15, 0))) |
    (p4_4 + theme(plot.margin = margin(-3, -25, -15, 0))) |
    p4_4_4
) + 
  plot_layout(guides = "collect") &
  theme(
    legend.position = "right",
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 8)
  )

# Display
combined_plot

# For better rendering when saving
ggsave("combined_heatmaps.png", 
       combined_plot,
       width = 13, 
       height = 10,
       dpi = 300,
       bg = "transparent")





  ##################################
 #        BOXPLOTS CI WIDTH       #
##################################


# Load necessary libraries
library(ggplot2)
library(reshape2)
library(dplyr)
library(patchwork)

load("./Run-all-10-no-concentrations-upper100/Adaptation_all_10_no_concentration-upper100.RData")
#load("./Run-3pg-2pg-10-no-concentrations-upper100/Adaptation-3pg-2pg-10-no-concentrations-upper100.RData")
#load("./Run-2pg-pep-10-no-concentrations-upper100/Adaptation-2pg-pep-10-no-concentrations-upper100.RData")

#load("./Run-all-10-concentrations-upper100/Adaptation_all_10_concentration-upper100.RData")
#load("./Run-3pg-2pg-10-concentrations-upper100/Adaptation_3pg_2pg_10_concentration-upper100.RData")

#load("./Run-all-10-concentrations-5-90-upper100/Adaptation_all_10_concentration-5-90-upper100.RData")
#load("./Run-3pg-2pg-10-concentrations-5-90-upper100/Adaptation_3pg_2pg_10_concentration-5-90-upper100.RData")

#load("./Run-all-10-irreversible-upper100/Adaptation_all-irreversible-upper100.RData")
#load("./Run-3pg-2pg-10-irreversible-upper100/Adaptation_3pg-2pg-irreversible-upper100.RData")
#load("./Run-2pg-pep-10-irreversible-upper100/Adaptation_2pg-pep-irreversible-upper100.RData")


#Gotta do a 3D matrix to store some info for plotting
# Define the strains and metrics
strains <- c("WT", "dcp12", "dcp12::cp12", "dcp12::cp12dCysC", "dcp12::cp12dCysN", "dcp12::cp12dCysNC")
metrics <- c("opt", "mean", "median", "sd", "ci_2.5", "ci_97.5")

# Initialize a 3D array to store results for all 10 runs
# Dimensions: strains x metrics x runs
conf_array <- array(NA, dim = c(length(strains), length(metrics), 10),
                    dimnames = list(strains, metrics, paste0("Run_", 1:10)))

# Loop through the 10 runs and fill the array
for (r in 1:10) {
  for (i in 1:length(strains)) {
    strain <- strains[i]
    
    # Extract the results for the current strain and run
    if (strain == "WT") {
      conf_array[i, , r] <- result_vars2[[r]]$resF_wt$sens$summary["v1", ]
    } else if (strain == "dcp12") {
      conf_array[i, , r] <- result_vars2[[r]]$resF_dcp12$sens$summary["v1", ]
    } else if (strain == "dcp12::cp12") {
      conf_array[i, , r] <- result_vars2[[r]]$resF_dcp12_cp12$sens$summary["v1", ]
    } else if (strain == "dcp12::cp12dCysC") {
      conf_array[i, , r] <- result_vars2[[r]]$resF_dcp12_cp12dCysC$sens$summary["v1", ]
    } else if (strain == "dcp12::cp12dCysN") {
      conf_array[i, , r] <- result_vars2[[r]]$resF_dcp12_cp12dCysN$sens$summary["v1", ]
    } else if (strain == "dcp12::cp12dCysNC") {
      conf_array[i, , r] <- result_vars2[[r]]$resF_dcp12_cp12dCysNC$sens$summary["v1", ]
    }
  }
}



# Extract data for plotting
plotci_data <- data.frame(
  Strain = character(),
  Run = numeric(),
  Median = numeric(),
  Mean = numeric(),
  CI_lower = numeric(),
  CI_upper = numeric()
)

# Loop through the 10 runs and strains to extract data
for (r in 1:10) {
  for (i in 1:length(strains)) {
    strain <- strains[i]
    
    # Extract the mean, ci_2.5, and ci_97.5 for the current strain and run
    if (strain == "WT") {
      mean_val <- result_vars2[[r]]$resF_wt$sens$summary["v1", "mean"]
      median_val <- result_vars2[[r]]$resF_wt$sens$summary["v1", "median"]
      ci_lower <- result_vars2[[r]]$resF_wt$sens$summary["v1", "ci_2.5"]
      ci_upper <- result_vars2[[r]]$resF_wt$sens$summary["v1", "ci_97.5"]
    } else if (strain == "dcp12") {
      mean_val <- result_vars2[[r]]$resF_dcp12$sens$summary["v1", "mean"]
      median_val <- result_vars2[[r]]$resF_dcp12$sens$summary["v1", "median"]
      ci_lower <- result_vars2[[r]]$resF_dcp12$sens$summary["v1", "ci_2.5"]
      ci_upper <- result_vars2[[r]]$resF_dcp12$sens$summary["v1", "ci_97.5"]
    } else if (strain == "dcp12::cp12") {
      mean_val <- result_vars2[[r]]$resF_dcp12_cp12$sens$summary["v1", "mean"]
      median_val <- result_vars2[[r]]$resF_dcp12_cp12$sens$summary["v1", "median"]
      ci_lower <- result_vars2[[r]]$resF_dcp12_cp12$sens$summary["v1", "ci_2.5"]
      ci_upper <- result_vars2[[r]]$resF_dcp12_cp12$sens$summary["v1", "ci_97.5"]
    } else if (strain == "dcp12::cp12dCysC") {
      mean_val <- result_vars2[[r]]$resF_dcp12_cp12dCysC$sens$summary["v1", "mean"]
      median_val <- result_vars2[[r]]$resF_dcp12_cp12dCysC$sens$summary["v1", "median"]
      ci_lower <- result_vars2[[r]]$resF_dcp12_cp12dCysC$sens$summary["v1", "ci_2.5"]
      ci_upper <- result_vars2[[r]]$resF_dcp12_cp12dCysC$sens$summary["v1", "ci_97.5"]
    } else if (strain == "dcp12::cp12dCysN") {
      mean_val <- result_vars2[[r]]$resF_dcp12_cp12dCysN$sens$summary["v1", "mean"]
      median_val <- result_vars2[[r]]$resF_dcp12_cp12dCysN$sens$summary["v1", "median"]
      ci_lower <- result_vars2[[r]]$resF_dcp12_cp12dCysN$sens$summary["v1", "ci_2.5"]
      ci_upper <- result_vars2[[r]]$resF_dcp12_cp12dCysN$sens$summary["v1", "ci_97.5"]
    } else if (strain == "dcp12::cp12dCysNC") {
      mean_val <- result_vars2[[r]]$resF_dcp12_cp12dCysNC$sens$summary["v1", "mean"]
      median_val <- result_vars2[[r]]$resF_dcp12_cp12dCysNC$sens$summary["v1", "median"]
      ci_lower <- result_vars2[[r]]$resF_dcp12_cp12dCysNC$sens$summary["v1", "ci_2.5"]
      ci_upper <- result_vars2[[r]]$resF_dcp12_cp12dCysNC$sens$summary["v1", "ci_97.5"]
    }
    
    # Add the data to the plot_data data frame
    plotci_data <- rbind(plotci_data, data.frame(
      Strain = strain,
      Run = r,
      Mean = mean_val,
      Median = median_val,
      CI_lower = ci_lower,
      CI_upper = ci_upper
    ))
  }
}


# Calculate CI width and add as new column
plotci_data$CI_width <- plotci_data$CI_upper - plotci_data$CI_lower

# Add yes/no column (50% threshold of Mean flux)
plotci_data$CI_acceptable <- ifelse(plotci_data$CI_width < 0.5 * plotci_data$Mean, "yes", "no")

# Coefficient of variation CV < 0.2 suggests high reliability
CV <- plotci_data %>%  
  group_by(Strain) %>%  
  summarise(CV = sd(Mean) / mean(Mean))

hits <- plotci_data %>%  
  group_by(Strain) %>%
  summarize(
    yes_count = sum(CI_acceptable == "yes"))

# View the modified dataframe
plotci_data

################### PAY ATTENTION HERE ##################
# To save later
plotci_data_4_4_4 <- plotci_data
CV_4_4_4 <- CV
hits_4_4_4 <- hits

# Define your custom color palette
col_factor <- c("WT" = "#1f77b4", 
                "dcp12" = "#ff7f0e", 
                "dcp12::cp12" = "#2ca02c", 
                "dcp12::cp12dCysC" = "#d62728", 
                "dcp12::cp12dCysN" = "#9467bd", 
                "dcp12::cp12dCysNC" = "#8c564b")

create_boxplot <- function(plotci_data, title, caption) {
  p <- ggplot(plotci_data, aes(x = factor(Run), y = Mean, color = Strain)) +
    geom_boxplot(
      aes(ymin = CI_lower, lower = CI_lower, middle = Median, 
          upper = CI_upper, ymax = CI_upper),
      stat = "identity",
      fill = "white",
      width = 0.6
    ) +
    # Add mean points (kept red for consistency)
    geom_point(
      aes(y = Mean),
      shape = 18,
      size = 4,
      color = "gold2"
    ) +
    # Add CI range lines (kept black for visibility)
    geom_errorbar(
      aes(ymin = CI_lower, ymax = CI_upper),
      width = 0.2,
      color = "black"
    ) +
    facet_wrap(~Strain, scales = "free_y", ncol = 3) +
    labs(
      title = title,
      x = "",
      y = "Flux with 95% CI range",
      caption = caption
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5),
      panel.spacing = unit(0.5, "lines"),
      legend.position = "bottom"
    ) +
    # Apply custom color scale
    scale_color_manual(values = col_factor)
  print(p)
}



# Create plots (only left column shows strain names)
p1 <- create_boxplot(plotci_data_1, "3PGA ↔ 2PGA ↔ PEP → ∅ No Conc 0-90 min.","")
p1_1 <- create_boxplot(plotci_data_1_1, "3PGA ↔ 2PGA → ∅ No Conc 0-90 min.","")
p1_1_1 <- create_boxplot(plotci_data_1_1_1, "2PGA ↔ PEP → ∅ No Conc 0-90 min.", "Boxes show 95% CI range (lower to upper)\nYellow diamond = Mean flux\nMiddle line = Median flux")

p2 <- create_boxplot(plotci_data_2, "3PGA ↔ 2PGA ↔ PEP → ∅ Conc 0-90 min.","")
p2_2 <- create_boxplot(plotci_data_2_2, "3PGA ↔ 2PGA → ∅ Conc 0-90 min.", "Boxes show 95% CI range (lower to upper)\nYellow diamond = Mean flux\nMiddle line = Median flux")

p3 <- create_boxplot(plotci_data_3, "3PGA ↔ 2PGA ↔ PEP → ∅ Conc 5-90 min.","")
p3_3 <- create_boxplot(plotci_data_3_3, "3PGA ↔ 2PGA → ∅ Conc 5-90 min.", "Boxes show 95% CI range (lower to upper)\nYellow diamond = Mean flux\nMiddle line = Median flux")

p4 <- create_boxplot(plotci_data_4, "3PGA → 2PGA → PEP → ∅ Conc 0-90 min.","")
p4_4 <- create_boxplot(plotci_data_4_4, "3PGA → 2PGA → ∅ Conc 0-90 min","")
p4_4_4 <- create_boxplot(plotci_data_4_4_4, "2PGA → PEP → ∅ Conc 0-90 min", "Boxes show 95% CI range (lower to upper)\nYellow diamond = Mean flux\nMiddle line = Median flux")




# No Conc 0-90

# Combine with a single caption
Scenario1 <- (p1 / p1_1 / p1_1_1) +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    plot.title = element_text(size = 10),
    plot.margin = margin(0,0,0,0)
  )
Scenario1


# Conc 0-90
Scenario2 <- (p2 / p2_2) +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    plot.title = element_text(size = 10),
    plot.margin = margin(0,0,0,0)
  )
Scenario2


# Conc 5-90
Scenario3 <- ( p3 / p3_3 ) +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    plot.title = element_text(size = 10),
    plot.margin = margin(0,0,0,0)
  )
Scenario3

# Conc irrev 0-90
Scenario4 <- ( p4 / p4_4 /  p4_4_4 ) +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    plot.title = element_text(size = 10),
    plot.margin = margin(0,0,0,0)
  )
Scenario4


# Keep these objects (won't be deleted)
keep_objects <- c("^plotci_data_", "^CV_", "^p[0-9]", "^hits", "keep_objects", "^Scenario[0-9]", "^case_", "^flags_", "merged")


# Remove objects NOT matching any pattern
rm(list = ls()[!grepl(paste(keep_objects, collapse = "|"), ls())])

save.image("NewBoxplots.RData")

merged_flux <- read_excel("./Merged_Max_Min_Flux.xlsx")


# Step 1: Calculate the true lowest min_flux across all metabolites
merged_flux <- merged_flux %>%
  rowwise() %>%  # Ensure calculations are done row-by-row
  mutate(
    lowest_min_flux = min(min_flux_3PGA, min_flux_2PGA, min_flux_PEP, na.rm = TRUE)
  ) %>%
  ungroup()

# Used the lowest of 2PGA as a flagger to detect false positives of reliable hits
# Step 2: Merge with plot data and flag
case_4_4_4 <- plotci_data_4_4_4 %>%
  left_join(
    merged_flux %>% select(Strain, lowest_min_flux),
    by = "Strain"
  ) %>%
  mutate(
    flag = ifelse(is.na(lowest_min_flux), FALSE, Mean < lowest_min_flux)
  )

# Step 3: Format labels with * for flagged cells
case_4_4_4 <- case_4_4_4 %>%
  mutate(
    Run = factor(Run, levels = 1:10),  # Force Run as ordered factor
    CI_formatted = case_when(
      is.na(Mean) ~ NA_character_,
      CI_acceptable == "yes" & flag ~ sprintf("<b>%.3f*</b>", CI_width),
      CI_acceptable == "yes" ~ sprintf("<b>%.3f</b>", CI_width),
      flag ~ sprintf("%.3f*", CI_width),
      TRUE ~ sprintf("%.3f", CI_width)
    )
  )

flags_1_1_1 <- case_1_1_1 %>%  
  group_by(Strain) %>%
  summarize(
    TRUE_count = sum(flag == "TRUE" & CI_acceptable == "yes"))

case_1 %>%  
  group_by(Strain) %>%
  summarize(
    TRUE_count = sum(flag == "TRUE" & CI_acceptable == "yes"))

flags_4_4_4
