# Adapted from:
#
# millard@insa-toulouse.fr
#
# IsoSim: stoichiometric, kinetic and isotopic modeling of metabolic systems
#
# https://github.com/MetaSys-LISBP/IsoSim
#
# https://doi.org/10.1101/735308
#
# Copyright 2019, INRA, France
# License: GNU General Public License v3 (see license.txt for details)

################### Run from 3PG <-> 2PG ###################

cat("
    This code calculates the flux through 3PG <-> 2PG in 6 Synechocystis strains
    \n \n")

####################################
cat("\n   ... Initialize R environment ...\n\n")
####################################

# Clean workspace
#rm(list= ls())

# Path of the prenyl pyrophosphate pathway model
#load("C:/Users/coron/Documents/Thesis/isosim/models/thesis_pathway/Run_3PG_2PG/3pg_2pg_fit.RData")
#library(Matrix)
#setwd("C:/Users/coron/Documents/Thesis/isosim/models/thesis_pathway/")
setwd("/home/coronamacias/Thesis/IsoSim/models/thesis_pathway/")

# Get current directory
wd <- getwd()

# Load IsoSim
setwd(file.path(dirname(dirname(wd)), "isosim"))
source("isosim.R")

# Go back to working directory and create "res" folder to store the results
setwd(wd)
if (!file.exists("Run-3pg-2pg-10-concentrations-upper100")){
  dir.create(file.path(wd, "Run-3pg-2pg-10-concentrations-upper100"))
}
setwd(file.path(wd, "Run-3pg-2pg-10-concentrations-upper100"))


#########################
### GLOBAL PARAMETERS ###
#########################

# number of cores to use in parallel, i.e. at most how many child processes will be run simultaneously
# single-core version can be used by setting 'numCores' to NULL
numCores <- 10

# number of iterations for Monte Carlo sensitivity analysis
niter <- 4

cat("\n Modified: concentrations are updated according to empirical data \n\n")
####################################
cat("\n   ... Construct the desired pathway (3PG <-> 2PG) ...\n\n")
####################################

rxn <- list(r1f = list("R"=c("PG2", "PG3"), "C"=c(1, -1),     "E"="v1",   "T"=c("ABC", "ABC")),
            r1r = list("R"=c("PG2", "PG3"), "C"=c(-1, 1),     "E"="v1",   "T"=c("ABC", "ABC")),
            r2 = list("R"=c("PG2"), "C"=c(-1),        "E"="v1",  "T"=c("ABC")))

net <- net2mat(rxn)


####################################
cat("\n   ... fit label inputs ...\n\n")
####################################

# measurement times (in min)
times <- c(0,5,10,15,30,60,90)
noise <- 0.0005
# 13C-enrichments of PG3 (label input) in strains "dcp12", "dcp12::cp12", "dcp12::cp12dCysC", "dcp12::cp12dCysN", "dcp12::cp12dCysNC", "WT K"
original_data_input <- cbind("WTK_PG3_1-2-3"=c(0, 0, 0, 0, 0, 0, 0),
                             "dcp12_PG3_1-2-3"=c(0,0,0.00397159,0.011314637,0.062717565,0.241552885,0.337020845),
                             "dcp12cp12_PG3_1-2-3"=c(0, 0, 0, 0, 0, 0.001792413, 0.00084248),
                             "dcp12cp12dCysC_PG3_1-2-3"=c(0, 0, 0, 0.00304893, 0.050659085, 0.0441527275, 0.0458769425),
                             "dcp12cp12dCysN_PG3_1-2-3"=c(0, 0, 0, 0, 0, 0, 0),
                             "dcp12cp12dCysNC_PG3_1-2-3"=c(0, 0, 0, 0.0161069625, 0.077143535, 0.219639593333333, 0.233304856666667))

data_input <- apply(original_data_input, 2, function(col) { #Add positive noise
  col + abs(rnorm(length(col), mean = 0, sd = noise))  })

# fit labeling dynamics of 3PG to define label input
enr_in <- fit_label_input(data_input, t=times, file="res_fit_enr", mc.cores=numCores)

# loop here check the fit, chi test, use noise within sd check against original 

# If there is time, do the presentation graphs smoother, point from 0 to 90, not just eval in the 7 measurements.


####################################
cat("\n   ... Calculate fluxes ...\n\n")
####################################

####################################
cat("\n      ... WT K strain\n\n")
####################################

# 13C-enrichments of metabolic intermediates (PG2)
enr_PG2_wt <- c(0, 0, 0, 0, 0, 0, 0)
data_exp_wt <- cbind(times, enr_PG2_wt)
colnames(data_exp_wt) <- c("time", "PG2_1-2-3")


# define analytical functions of label input
# here we use functions from the best fit
anFun_fit_wt <- list("PG3_1-2-3-M0" = enr_in$anFun[[1]], #Changed ipp and adjusted to what there is in enr_in$anFun ★
                     "PG3_1-2-3-M1" = enr_in$anFun[[2]],
                     "PG3_1-2-3-M2" = enr_in$anFun[[3]],
                     "PG3_1-2-3-M3" = enr_in$anFun[[4]])

# initial parameters
#    metabolite concentrations
meta_conc_wt <- c("PG3"= 1231.698242, "PG2"=5.479722338)
#    fluxes
kp_wt <- c("v1"=0)
#    sd of measurements
sd_meas_wt <- list(iso=0.5, conc=c("PG2"=0.023380))
#    names of parameters (concentrations & fluxes) to estimate  
te_wt <- c("v1", "PG2")  # Here is somehow less parameters in the original example
#    upper & lower constraints on parameters
te_upc_wt <- c("v1"=10, "PG2"=1)
te_loc_wt <- c("v1"=1e-4, "PG2"=1e-4)

# run flux calculation
resF_wt <- fit(net        = net,
               times      = times,
               kp         = kp_wt,
               to_est     = te_wt,
               te_upc     = te_upc_wt[te_wt],
               te_loc     = te_loc_wt[te_wt],
               te_upc_det = NULL,
               te_loc_det = NULL,
               eq_det     = NULL,
               data_meas  = list(iso=data_exp_wt, conc=c("PG2"=0.023380)),
               sd_meas    = sd_meas_wt,
               meta_conc  = meta_conc_wt,
               iso_conc   = NULL,
               mode       = "enrichments",
               anFun      = anFun_fit_wt,
               trf        = NULL,
               events     = NULL,
               p          = 0.0,
               subDir     = "fit_wt",
               nlsic_ctrl = list(errx=1.e-6, maxit=50, trace=1, ci=list(report=T), history=F, adaptbt=T),
               method     = "FORTRAN",
               niter      = niter)

# display estimated parameters
resF_wt$sens$summary


####################################
cat("\n      ... dcp12 strain ...     \n\n")
####################################

# 13C-enrichments of metabolic intermediates (PG2A, PEP)
enr_PG2_dcp12 <- c(0, 0, 0.0010236875, 0.005175985, 0.1018379825, 0.2487054925, 0.3470264875)
data_exp_dcp12 <- cbind(times, enr_PG2_dcp12)
colnames(data_exp_dcp12) <- c("time", "PG2_1-2-3")


# define analytical functions of label input
# here we use functions from the best fit
anFun_fit_dcp12 <- list("PG3_1-2-3-M0" = enr_in$anFun[[5]], #Changed ipp and adjusted to what there is in enr_in$anFun ★
                        "PG3_1-2-3-M1" = enr_in$anFun[[6]],
                        "PG3_1-2-3-M2" = enr_in$anFun[[7]],
                        "PG3_1-2-3-M3" = enr_in$anFun[[8]])

# initial parameters
#    metabolite concentrations
meta_conc_dcp12 <- c("PG3"=1290.13851733636, "PG2"=8.23670666621142)
#    fluxes
kp_dcp12 <- c("v1"=1)
#    sd of measurements
sd_meas_dcp12 <- list(iso=sd(enr_PG2_dcp12), conc=c("PG2"=0.022175))
#    names of parameters (concentrations & fluxes) to estimate  
te_dcp12 <- c("v1", "PG2")  # Here is somehow less parameters in the original example
#    upper & lower constraints on parameters
te_upc_dcp12 <- c("v1"=100, "PG2"=1)
te_loc_dcp12 <- c("v1"=1e-4, "PG2"=1e-4)

# run flux calculation
resF_dcp12 <- fit(net        = net,
                  times      = times,
                  kp         = kp_dcp12,
                  to_est     = te_dcp12,
                  te_upc     = te_upc_dcp12[te_dcp12],
                  te_loc     = te_loc_dcp12[te_dcp12],
                  te_upc_det = NULL,
                  te_loc_det = NULL,
                  eq_det     = NULL,
                  data_meas  = list(iso=data_exp_dcp12, conc=c("PG2"=0.022175)),
                  sd_meas    = sd_meas_dcp12,
                  meta_conc  = meta_conc_dcp12,
                  iso_conc   = NULL,
                  mode       = "enrichments",
                  anFun      = anFun_fit_dcp12,
                  trf        = NULL,
                  events     = NULL,
                  p          = 0.0,
                  subDir     = "fit_dcp12",
                  nlsic_ctrl = list(errx=1.e-6, maxit=50, trace=1, ci=list(report=T), history=F, adaptbt=T),
                  method     = "FORTRAN",
                  niter      = niter)

# display estimated parameters
resF_dcp12$sens$summary


####################################
cat("\n      ... dcp12::cp12 strain ... \n\n")
####################################

# 13C-enrichments of metabolic intermediates (PG2A, PEP)
enr_PG2_dcp12_cp12 <- c(0, 0, 0, 0, 0, 0.000253238, 0.003050705)
data_exp_dcp12_cp12 <- cbind(times, enr_PG2_dcp12_cp12)
colnames(data_exp_dcp12_cp12) <- c("time", "PG2_1-2-3")


# define analytical functions of label input
# here we use functions from the best fit
anFun_fit_dcp12_cp12 <- list("PG3_1-2-3-M0" = enr_in$anFun[[9]], #Changed ipp and adjusted to what there is in enr_in$anFun ★
                             "PG3_1-2-3-M1" = enr_in$anFun[[10]],
                             "PG3_1-2-3-M2" = enr_in$anFun[[11]],
                             "PG3_1-2-3-M3" = enr_in$anFun[[12]])

# initial parameters
#    metabolite concentrations
meta_conc_dcp12_cp12 <- c("PG3"=1045.01809443301, "PG2"=5.74411473678459) 
#    fluxes
kp_dcp12_cp12 <- c("v1"=0)
#    sd of measurements
sd_meas_dcp12_cp12 <- list(iso=sd(enr_PG2_dcp12_cp12), conc=c("PG2"=0.013020))
#    names of parameters (concentrations & fluxes) to estimate  
te_dcp12_cp12 <- c("v1", "PG2")  # Here is somehow less parameters in the original example
#    upper & lower constraints on parameters
te_upc_dcp12_cp12 <- c("v1"=10, "PG2"=1) #changed upper limit. Before got error by just changing kp_dcp12_cp12 to 0
te_loc_dcp12_cp12 <- c("v1"=1e-4, "PG2"=1e-4)

# run flux calculation
resF_dcp12_cp12 <- fit(net        = net,
                       times      = times,
                       kp         = kp_dcp12_cp12,
                       to_est     = te_dcp12_cp12,
                       te_upc     = te_upc_dcp12_cp12[te_dcp12_cp12],
                       te_loc     = te_loc_dcp12_cp12[te_dcp12_cp12],
                       te_upc_det = NULL,
                       te_loc_det = NULL,
                       eq_det     = NULL,
                       data_meas  = list(iso=data_exp_dcp12_cp12, conc=c("PG2"=0.013020)),
                       sd_meas    = sd_meas_dcp12_cp12,
                       meta_conc  = meta_conc_dcp12_cp12,
                       iso_conc   = NULL,
                       mode       = "enrichments",
                       anFun      = anFun_fit_dcp12_cp12,
                       trf        = NULL,
                       events     = NULL,
                       p          = 0.0,
                       subDir     = "fit_dcp12_cp12",
                       nlsic_ctrl = list(errx=1.e-6, maxit=50, trace=1, ci=list(report=T), history=F, adaptbt=T),
                       method     = "FORTRAN",
                       niter      = niter)

# display estimated parameters
resF_dcp12_cp12$sens$summary


####################################
cat("\n      ... dcp12::cp12dCysC strain\n\n")
####################################

# 13C-enrichments of metabolic intermediates (PG2A, PEP)
enr_PG2_dcp12_cp12dCysC <- c(0, 0, 0, 0.0067180525, 0.0495900833333333, 0.0565733, 0.0485604875)
data_exp_dcp12_cp12dCysC <- cbind(times, enr_PG2_dcp12_cp12dCysC)
colnames(data_exp_dcp12_cp12dCysC) <- c("time", "PG2_1-2-3")


# define analytical functions of label input
# here we use functions from the best fit
anFun_fit_dcp12_cp12dCysC <- list("PG3_1-2-3-M0" = enr_in$anFun[[13]], #Changed ipp and adjusted to what there is in enr_in$anFun ★
                                  "PG3_1-2-3-M1" = enr_in$anFun[[14]],
                                  "PG3_1-2-3-M2" = enr_in$anFun[[15]],
                                  "PG3_1-2-3-M3" = enr_in$anFun[[16]])

# initial parameters
#    metabolite concentrations
meta_conc_dcp12_cp12dCysC <- c("PG3"=859.155016461378, "PG2"=10.611146069799) 
#    fluxes
kp_dcp12_cp12dCysC <- c("v1"=1)
#    sd of measurements
sd_meas_dcp12_cp12dCysC <- list(iso=sd(enr_PG2_dcp12_cp12dCysC), conc=c("PG2"=0.017141))
#    names of parameters (concentrations & fluxes) to estimate  
te_dcp12_cp12dCysC <- c("v1", "PG2")  # Here is somehow less parameters in the original example
#    upper & lower constraints on parameters
te_upc_dcp12_cp12dCysC <- c("v1"=10, "PG2"=1)
te_loc_dcp12_cp12dCysC <- c("v1"=1e-4, "PG2"=1e-4)

# run flux calculation
resF_dcp12_cp12dCysC <- fit(net        = net,
                            times      = times,
                            kp         = kp_dcp12_cp12dCysC,
                            to_est     = te_dcp12_cp12dCysC,
                            te_upc     = te_upc_dcp12_cp12dCysC[te_dcp12_cp12dCysC],
                            te_loc     = te_loc_dcp12_cp12dCysC[te_dcp12_cp12dCysC],
                            te_upc_det = NULL,
                            te_loc_det = NULL,
                            eq_det     = NULL,
                            data_meas  = list(iso=data_exp_dcp12_cp12dCysC, conc=c("PG2"=0.017141)),
                            sd_meas    = sd_meas_dcp12_cp12dCysC,
                            meta_conc  = meta_conc_dcp12_cp12dCysC,
                            iso_conc   = NULL,
                            mode       = "enrichments",
                            anFun      = anFun_fit_dcp12_cp12dCysC,
                            trf        = NULL,
                            events     = NULL,
                            p          = 0.0,
                            subDir     = "fit_dcp12_cp12dCysC",
                            nlsic_ctrl = list(errx=1.e-6, maxit=50, trace=1, ci=list(report=T), history=F, adaptbt=T),
                            method     = "FORTRAN",
                            niter      = niter)

# display estimated parameters
resF_dcp12_cp12dCysC$sens$summary


####################################
cat("\n      ... dcp12::cp12dCysN strain\n\n")
####################################

# 13C-enrichments of metabolic intermediates (PG2A, PEP)
enr_PG2_dcp12_cp12dCysN <- c(0, 0, 0, 0, 0, 0, 0)
data_exp_dcp12_cp12dCysN <- cbind(times, enr_PG2_dcp12_cp12dCysN)
colnames(data_exp_dcp12_cp12dCysN) <- c("time", "PG2_1-2-3")


# define analytical functions of label input
# here we use functions from the best fit
anFun_fit_dcp12_cp12dCysN <- list("PG3_1-2-3-M0" = enr_in$anFun[[17]], #Changed ipp and adjusted to what there is in enr_in$anFun ★
                                  "PG3_1-2-3-M1" = enr_in$anFun[[18]],
                                  "PG3_1-2-3-M2" = enr_in$anFun[[19]],
                                  "PG3_1-2-3-M3" = enr_in$anFun[[20]])

# initial parameters
#    metabolite concentrations
meta_conc_dcp12_cp12dCysN <- c("PG3"=966.898635494004, "PG2"=6.03444476511343) 
#    fluxes
kp_dcp12_cp12dCysN <- c("v1"=1)
#    sd of measurements
sd_meas_dcp12_cp12dCysN <- list(iso=0.5, conc=c("PG2"=0.013923))
#    names of parameters (concentrations & fluxes) to estimate  
te_dcp12_cp12dCysN <- c("v1", "PG2")  # Here is somehow less parameters in the original example
#    upper & lower constraints on parameters
te_upc_dcp12_cp12dCysN <- c("v1"=10, "PG2"=1)
te_loc_dcp12_cp12dCysN <- c("v1"=1e-4, "PG2"=1e-4)

# run flux calculation
resF_dcp12_cp12dCysN <- fit(net        = net,
                            times      = times,
                            kp         = kp_dcp12_cp12dCysN,
                            to_est     = te_dcp12_cp12dCysN,
                            te_upc     = te_upc_dcp12_cp12dCysN[te_dcp12_cp12dCysN],
                            te_loc     = te_loc_dcp12_cp12dCysN[te_dcp12_cp12dCysN],
                            te_upc_det = NULL,
                            te_loc_det = NULL,
                            eq_det     = NULL,
                            data_meas  = list(iso=data_exp_dcp12_cp12dCysN, conc=c("PG2"=0.013923)),
                            sd_meas    = sd_meas_dcp12_cp12dCysN,
                            meta_conc  = meta_conc_dcp12_cp12dCysN,
                            iso_conc   = NULL,
                            mode       = "enrichments",
                            anFun      = anFun_fit_dcp12_cp12dCysN,
                            trf        = NULL,
                            events     = NULL,
                            p          = 0.0,
                            subDir     = "fit_dcp12_cp12dCysN",
                            nlsic_ctrl = list(errx=1.e-6, maxit=50, trace=1, ci=list(report=T), history=F, adaptbt=T),
                            method     = "FORTRAN",
                            niter      = niter)

# display estimated parameters
resF_dcp12_cp12dCysN$sens$summary


####################################
cat("\n      ... dcp12::cp12dCysNC strain\n\n")
####################################

# 13C-enrichments of metabolic intermediates (PG2A, PEP)
enr_PG2_dcp12_cp12dCysNC <- c(0, 0, 0, 0.011180058, 0.076752798, 0.204328553, 0.24085604)
data_exp_dcp12_cp12dCysNC <- cbind(times, enr_PG2_dcp12_cp12dCysNC)
colnames(data_exp_dcp12_cp12dCysNC) <- c("time", "PG2_1-2-3")


# define analytical functions of label input
# here we use functions from the best fit
anFun_fit_dcp12_cp12dCysNC <- list("PG3_1-2-3-M0" = enr_in$anFun[[21]], #Changed ipp and adjusted to what there is in enr_in$anFun ★
                                   "PG3_1-2-3-M1" = enr_in$anFun[[22]],
                                   "PG3_1-2-3-M2" = enr_in$anFun[[23]],
                                   "PG3_1-2-3-M3" = enr_in$anFun[[24]])

# initial parameters
#    metabolite concentrations
meta_conc_dcp12_cp12dCysNC <- c("PG3"=954.676194927616, "PG2"=5.50614272562167) 
#    fluxes
kp_dcp12_cp12dCysNC <- c("v1"=1)
#    sd of measurements
sd_meas_dcp12_cp12dCysNC <- list(iso=sd(enr_PG2_dcp12_cp12dCysNC), conc=c("PG2"=0.015778))
#    names of parameters (concentrations & fluxes) to estimate  
te_dcp12_cp12dCysNC<- c("v1", "PG2")  # Here is somehow less parameters in the original example
#    upper & lower constraints on parameters
te_upc_dcp12_cp12dCysNC <- c("v1"=10, "PG2"=1)
te_loc_dcp12_cp12dCysNC <- c("v1"=1e-4, "PG2"=1e-4)

# run flux calculation
resF_dcp12_cp12dCysNC <- fit(net        = net,
                             times      = times,
                             kp         = kp_dcp12_cp12dCysNC,
                             to_est     = te_dcp12_cp12dCysNC,
                             te_upc     = te_upc_dcp12_cp12dCysNC[te_dcp12_cp12dCysNC],
                             te_loc     = te_loc_dcp12_cp12dCysNC[te_dcp12_cp12dCysNC],
                             te_upc_det = NULL,
                             te_loc_det = NULL,
                             eq_det     = NULL,
                             data_meas  = list(iso=data_exp_dcp12_cp12dCysNC, conc=c("PG2"=0.015778)),
                             sd_meas    = sd_meas_dcp12_cp12dCysNC,
                             meta_conc  = meta_conc_dcp12_cp12dCysNC,
                             iso_conc   = NULL,
                             mode       = "enrichments",
                             anFun      = anFun_fit_dcp12_cp12dCysNC,
                             trf        = NULL,
                             events     = NULL,
                             p          = 0.0,
                             subDir     = "fit_dcp12_cp12dCysNC",
                             nlsic_ctrl = list(errx=1.e-6, maxit=50, trace=1, ci=list(report=T), history=F, adaptbt=T),
                             method     = "FORTRAN",
                             niter      = niter)

# display estimated parameters
resF_dcp12_cp12dCysNC$sens$summary


####################################
cat("\n   ... plot results ...\n\n")
####################################

# colors for "WT", "dcp12", "dcp12::cp12", "dcp12::cp12dCysC", "dcp12::cp12dCysN and "dcp12::cp12dCysNC" datasets
col_f <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b")


# Small offset to avoid plotting issues with zero-only values
epsilon <- 1e-6

# Standard deviation for WT
sd_wt <- list("3PG"=c(0, 0, 0, 0, 0, 0, 0),
              "2PG"=c(0, 0, 0, 0, 0, 0, 0))


#################################################################################

# dcp12


# Standard deviation for dcp12
sd_dcp12 <- list(
  "3PG"=c(0, 0, 0.00794318, 0.010755359, 0.043902583, 0.03444038, 0.044244396),
  "2PG"=c(0, 0, 0.002047375, 0.007319948, 0.034679842, 0.043441443, 0.075407513)
)


################################################################################################
# dcp1::cp12


# Standard deviation for dcp12::cp12
sd_dcp12_cp12 <- list(
  "3PG"=c(0, 0, 0, 0, 0, 0.003584825,	0.00168496),
  "2PG"=c(0, 0,	0, 0,	0, 0.000506475,	0.00610141)
)


##############################################################################################
# dcp12::cp12dCysC


# Standard deviation for dcp12::cp12
sd_dcp12_cp12dCysC <- list(
  "3PG"=c(0, 0, 0, 0.003532703, 0.009230331, 0.030750634, 0.009451074),
  "2PG"=c(0, 0, 0, 0.007980192, 0.01988221, 0.011417744, 0.013267179)
)


##############################################################################################
# dcp12::cp12dCysN


# Standard deviation for dcp12::cp12
sd_dcp12_cp12dCysN <- list(
  "3PG"=c(0, 0, 0, 0, 0, 0, 0),
  "2PG"=c(0, 0, 0, 0, 0, 0, 0)
)

##############################################################################################
# dcp12::cp12dCysNC


# Standard deviation for dcp12::cp12
sd_dcp12_cp12dCysNC <- list(
  "3PG"=c(0, 0, 0, 0.018721401, 0.030749119, 0.02480931, 0.050393998),
  "2PG"=c(0, 0, 0, 0.019012441, 0.029334965, 0.023493797, 0.037075651)
)

################################################################################
################################################################################

flux <- c(resF_wt$result$par["v1"],
          resF_dcp12$result$par["v1"],
          resF_dcp12_cp12$result$par["v1"],
          resF_dcp12_cp12dCysC$result$par["v1"],
          resF_dcp12_cp12dCysN$result$par["v1"],
          resF_dcp12_cp12dCysNC$result$par["v1"]
)
sd_flux<- c(resF_wt$sens$summary["v1", "sd"],
            resF_dcp12$sens$summary["v1", "sd"],
            resF_dcp12_cp12$sens$summary["v1", "sd"],
            resF_dcp12_cp12dCysC$sens$summary["v1", "sd"],
            resF_dcp12_cp12dCysN$sens$summary["v1", "sd"],
            resF_dcp12_cp12dCysNC$sens$summary["v1", "sd"])


turnover_2PG <- c(resF_wt$result$par["PG2"]/resF_wt$result$par["v1"],
                  resF_dcp12$result$par["PG2"]/resF_dcp12$result$par["v1"],
                  resF_dcp12_cp12$result$par["PG2"]/resF_dcp12_cp12$result$par["v1"],
                  resF_dcp12_cp12dCysC$result$par["PG2"]/resF_dcp12_cp12dCysC$result$par["v1"],
                  resF_dcp12_cp12dCysN$result$par["PG2"]/resF_dcp12_cp12dCysN$result$par["v1"],
                  resF_dcp12_cp12dCysNC$result$par["PG2"]/resF_dcp12_cp12dCysNC$result$par["v1"])

sd_turnover_2PG <- c(sqrt(resF_wt$sens$summary["v1", "sd"]**2+resF_wt$sens$summary["PG2", "sd"]**2),
                     sqrt(resF_dcp12$sens$summary["v1", "sd"]**2+resF_dcp12$sens$summary["PG2", "sd"]**2),
                     sqrt(resF_dcp12_cp12$sens$summary["v1", "sd"]**2+resF_dcp12_cp12$sens$summary["PG2", "sd"]**2),
                     sqrt(resF_dcp12_cp12dCysC$sens$summary["v1", "sd"]**2+resF_dcp12_cp12dCysC$sens$summary["PG2", "sd"]**2),
                     sqrt(resF_dcp12_cp12dCysN$sens$summary["v1", "sd"]**2+resF_dcp12_cp12dCysN$sens$summary["PG2", "sd"]**2),
                     sqrt(resF_dcp12_cp12dCysNC$sens$summary["v1", "sd"]**2+resF_dcp12_cp12dCysNC$sens$summary["PG2", "sd"]**2))



