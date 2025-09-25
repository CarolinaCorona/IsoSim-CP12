######################### Simulate 10 runs ################################

# Change accordingly to which scenario and system is being run
rm(list = ls())

old_time <- Sys.time()

# To allocate results
results_flux <- matrix(NA, nrow = 10, ncol = 6)
colnames(results_flux) <- c("WT", "dcp12", "dcp12::cp12", "dcp12::cp12dCysC",
                            "dcp12::cp12dCysN", "dcp12::cp12dCysNC")

results_sd_fluxes <- matrix(NA, nrow = 10, ncol = 6)
colnames(results_sd_fluxes) <- c("WT", "dcp12", "dcp12::cp12", "dcp12::cp12dCysC",
                            "dcp12::cp12dCysN", "dcp12::cp12dCysNC")

# In cases of partial system 3PGA<->2PGA, here we change PEP turnover for 2PGA turnover
results_PEP_turnover <- matrix(NA, nrow = 10, ncol = 6)
colnames(results_PEP_turnover) <- c("WT", "dcp12", "dcp12::cp12", "dcp12::cp12dCysC",
                                 "dcp12::cp12dCysN", "dcp12::cp12dCysNC")

results_sd_PEP_turnover <- matrix(NA, nrow = 10, ncol = 6)
colnames(results_sd_PEP_turnover) <- c("WT", "dcp12", "dcp12::cp12", "dcp12::cp12dCysC",
                                 "dcp12::cp12dCysN", "dcp12::cp12dCysNC")
# To allocate fittings
data_input_used <- list()

results_fitting <- list()

plotting_vars <- list()

result_vars <- list()

result_vars2 <- list()

# Change the NOTE accordingly to keeo better track of results
# Run the loop 10 times
for (RUN in 1:10) {
  cat("NOTE: This is the code with the concentrations updated and upper 100 bound for cp12 for 2pga pep metabolites from 0 to 90 mins with irreversible reactions!")
  cat("\n\n")
  print(RUN)
  cat("\n\n")
  print(results_flux)
  
  # Source your script to generate flux and results
  source(".IsoSim-CP12/Scenario 4/Adaptation IsoSim Concentration Irreversible 2PGA PEP.R") # Double check name of file
  
  # Store the vector in the i-th row of the matrix
  results_flux[RUN, ] <- flux
  
  # Store an_Fun fittings
  results_fitting[[RUN]] <- list(
    WT = anFun_fit_wt,
    dcp12 = anFun_fit_dcp12,
    dcp12_cp12 = anFun_fit_dcp12_cp12,
    dcp12_cp12dCysC = anFun_fit_dcp12_cp12dCysC,
    dcp12_cp12dCysN = anFun_fit_dcp12_cp12dCysN,
    dcp12_cp12dCysNC = anFun_fit_dcp12_cp12dCysNC
  )
  
  # Define a list of result variables
  result_vars <- list(
    resF_wt = resF_wt,
    resF_dcp12 = resF_dcp12,
    resF_dcp12_cp12 = resF_dcp12_cp12,
    resF_dcp12_cp12dCysC = resF_dcp12_cp12dCysC,
    resF_dcp12_cp12dCysN = resF_dcp12_cp12dCysN,
    resF_dcp12_cp12dCysNC = resF_dcp12_cp12dCysNC
  )
  
  result_vars2[[RUN]] <- list(
    resF_wt = resF_wt,
    resF_dcp12 = resF_dcp12,
    resF_dcp12_cp12 = resF_dcp12_cp12,
    resF_dcp12_cp12dCysC = resF_dcp12_cp12dCysC,
    resF_dcp12_cp12dCysN = resF_dcp12_cp12dCysN,
    resF_dcp12_cp12dCysNC = resF_dcp12_cp12dCysNC
  )
  
  # Store other variables for plotting
  plotting_vars[[RUN]] <- list(
    resF_wt = resF_wt$result$retres$sim,
    resF_dcp12 = resF_dcp12$result$retres$sim,
    resF_dcp12_cp12 = resF_dcp12_cp12$result$retres$sim,
    resF_dcp12_cp12dCysC = resF_dcp12_cp12dCysC$result$retres$sim,
    resF_dcp12_cp12dCysN = resF_dcp12_cp12dCysN$result$retres$sim,
    resF_dcp12_cp12dCysNC = resF_dcp12_cp12dCysNC$result$retres$sim
  )
  
  results_sd_fluxes[RUN, ] <- sd_flux
  
  results_PEP_turnover[RUN, ] <- turnover_PEP
  
  results_sd_PEP_turnover[RUN, ] <- sd_turnover_PEP
  
  data_input_used[[RUN]] <- list(data_input)
  
  # Add each result as a new sheet in the workbook
  for (name in names(result_vars)) {
    # Extract the variable
    data_to_write <- result_vars2[[RUN]][name]
    # Text file
    # Create a file name for each iteration
    file_name <- paste0("Full_results_", RUN, "_", name,".txt")
    
    # Save the list to a text file
    capture.output(print(data_to_write), file = file_name)
    
    # You can also append additional information to the file if needed
    cat("Results for iteration", RUN,"for ", name, "saved to file:", file_name, "\n")
  }
  
  # Print a message to indicate completion of the iteration
  print(paste0("Iteration number ", RUN, " DONE!"))
  
  # Save workspace
  save.image("Adaptation_2pg-pep-irreversible-upper100.RData") # Change if necessary
  
  #dev.off()
  
  setwd(wd)
  
}
################################################
cat("\n ***** Runs concluded :) ***** \n")
###############################################

cat("\n\n READY FOR PLOTTING \n\n")

setwd("/home/coronamacias/Thesis/IsoSim/models/thesis_pathway/Run-2pg-pep-10-irreversible-upper100") # Change accordingly

save.image("Adaptation_2pg-pep-irreversible-upper100.RData") # Change accordingly

new_time <- Sys.time()

cat("\n\n Checking time taken for this process: \n")
print(new_time - old_time)

print(results_flux)

setwd(wd)

