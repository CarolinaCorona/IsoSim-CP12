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

################### Run from 2PG <-> PEP ###################

cat("
    This code calculates the flux through 2PG -> PEP -> 0 in 6 Synechocystis strains
    \n \n")

####################################
cat("\n   ... Initialize R environment ...\n\n")
####################################


# Path of the prenyl pyrophosphate pathway model
#setwd("C:/Users/coron/Documents/Thesis/isosim/models/thesis_pathway/")
setwd("/home/coronamacias/Thesis/IsoSim/models/thesis_pathway/")

# Get current directory
wd <- getwd()

# Load IsoSim
setwd(file.path(dirname(dirname(wd)), "isosim"))
source("isosim.R")

# Go back to working directory and create "res" folder to store the results
setwd(wd)
if (!file.exists("Run-2pg-pep-10-irreversible-upper100")){
  dir.create(file.path(wd, "Run-2pg-pep-10-irreversible-upper100"))
}
setwd(file.path(wd, "Run-2pg-pep-10-irreversible-upper100"))


#########################
### GLOBAL PARAMETERS ###
#########################

# number of cores to use in parallel, i.e. at most how many child processes will be run simultaneously
# single-core version can be used by setting 'numCores' to NULL
numCores <- 12

# number of iterations for Monte Carlo sensitivity analysis
niter <- 4


####################################
cat("\n   ... Construct the desired pathway (2PG -> PEP -> 0) ...\n\n")
####################################


rxn <- list(r1 = list("R"=c("PEP", "PG2"),    "C"=c(1, -1),     "E"="v1",   "T"=c("ABC", "ABC")),
            r2 = list("R"=c("PEP"), "C"=c(-1),        "E"="v1",  "T"=c("ABC")))

net <- net2mat(rxn)


####################################
cat("\n   ... fit label inputs ...\n\n")
####################################

# measurement times (in min)
times <- c(0,5,10,15,30,60,90)

noise <- 0.0005

# 13C-enrichments of PG2 (label input) in strains "dcp12", "dcp12::cp12", "dcp12::cp12dCysC", "dcp12::cp12dCysN", "dcp12::cp12dCysNC", "WT K"
original_data_input <- cbind("WTK_PG2_1-2-3"=c(0, 0, 0, 0, 0, 0, 0),
                             "dcp12_PG2_1-2-3"=c(0, 0, 0.0010236875, 0.005175985, 0.1018379825, 0.2487054925, 0.3470264875),
                             "dcp12cp12_PG2_1-2-3"=c(0, 0, 0, 0, 0, 0.0002532375, 0.003050705),
                             "dcp12cp12dCysC_PG2_1-2-3"=c(0, 0, 0, 0.0067180525, 0.0495900833333333, 0.0565733, 0.0485604875),
                             "dcp12cp12dCysN_PG2_1-2-3"=c(0, 0, 0, 0, 0, 0, 0),
                             "dcp12cp12dCysNC_PG2_1-2-3"=c(0, 0, 0, 0.0111800575, 0.0767527975, 0.204328553333333, 0.24085604))

data_input <- apply(original_data_input, 2, function(col) { #Add positive noise
  col + abs(rnorm(length(col), mean = 0, sd = noise))  })

# fit labeling dynamics of IPP to define label input
enr_in <- fit_label_input(data_input, t=times, file="res_fit_enr", mc.cores=numCores)


####################################
cat("\n   ... calculate fluxes ...\n\n")
####################################

####################################
cat("\n      ... WT K strain\n\n")
####################################

# 13C-enrichments of metabolic intermediates (PEP)
enr_PEP_wt <- c(0, 0, 0, 0, 0, 0, 0)
data_exp_wt <- cbind(times, enr_PEP_wt)
colnames(data_exp_wt) <- c("time", "PEP_1-2-3")


# define analytical functions of label input
# here we use functions from the best fit
anFun_fit_wt <- list("PG2_1-2-3-M0" = enr_in$anFun[[1]], #Changed ipp and adjusted to what there is in enr_in$anFun ★
                     "PG2_1-2-3-M1" = enr_in$anFun[[2]],
                     "PG2_1-2-3-M2" = enr_in$anFun[[3]],
                     "PG2_1-2-3-M3" = enr_in$anFun[[4]])

# initial parameters
#    metabolite concentrations
meta_conc_wt <- c("PG2"=5.479722338, "PEP"=743.7367661) #Most os Sebastian initial met concentr. 1.100000e-01 or 1.100000e-02
#    fluxes
kp_wt <- c("v1"=0)
#    sd of measurements
sd_meas_wt <- list(iso=0.5, conc=c("PEP"=1.399425))
#    names of parameters (concentrations & fluxes) to estimate  
te_wt <- c("v1", "PEP")  # Here is somehow less parameters in the original example
#    upper & lower constraints on parameters
te_upc_wt <- c("v1"=10, "PEP"=10)
te_loc_wt <- c("v1"=1e-4, "PEP"=1e-1)

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
               data_meas  = list(iso=data_exp_wt, conc=c("PEP"=1.399425)),
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
cat("\n      ... dcp12 strain\n\n")
####################################

# 13C-enrichments of metabolic intermediates (PEP)
enr_PEP_dcp12 <- c(0, 0, 0.007159935, 0.0217724333333333, 0.10046094, 0.2625095, 0.343587925)
data_exp_dcp12 <- cbind(times, enr_PEP_dcp12)
colnames(data_exp_dcp12) <- c("time", "PEP_1-2-3")


# define analytical functions of label input
# here we use functions from the best fit
anFun_fit_dcp12 <- list("PG2_1-2-3-M0" = enr_in$anFun[[5]], #Changed ipp and adjusted to what there is in enr_in$anFun ★
                        "PG2_1-2-3-M1" = enr_in$anFun[[6]],
                        "PG2_1-2-3-M2" = enr_in$anFun[[7]],
                        "PG2_1-2-3-M3" = enr_in$anFun[[8]])

# initial parameters
#    metabolite concentrations
meta_conc_dcp12 <- c("PG2"=8.23670666621142, "PEP"=870.883044294975) #Most os Sebastian initial met concentr. 1.100000e-01 or 1.100000e-02
#    fluxes
kp_dcp12 <- c("v1"=1)
#    sd of measurements
sd_meas_dcp12 <- list(iso=sd(enr_PEP_dcp12), conc=c("PEP"=1.750731279))
#    names of parameters (concentrations & fluxes) to estimate  
te_dcp12 <- c("v1", "PEP")  # Here is somehow less parameters in the original example
#    upper & lower constraints on parameters
te_upc_dcp12 <- c("v1"=100, "PEP"=15)
te_loc_dcp12 <- c("v1"=1e-4, "PEP"=1e-1)

# run flux calculation
resF_dcp12 <- fit(net     = net,
                  times      = times,
                  kp         = kp_dcp12,
                  to_est     = te_dcp12,
                  te_upc     = te_upc_dcp12[te_dcp12],
                  te_loc     = te_loc_dcp12[te_dcp12],
                  te_upc_det = NULL,
                  te_loc_det = NULL,
                  eq_det     = NULL,
                  data_meas  = list(iso=data_exp_dcp12, conc=c("PEP"=1.750731279)),
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
cat("\n      ... dcp12::cp12 strain\n\n")
####################################

# 13C-enrichments of metabolic intermediates (PEP)
enr_PEP_dcp12_cp12 <- c(0, 0, 0, 0, 0, 0.004714665, 0.0074160425)
data_exp_dcp12_cp12 <- cbind(times, enr_PEP_dcp12_cp12)
colnames(data_exp_dcp12_cp12) <- c("time", "PEP_1-2-3")


# define analytical functions of label input
# here we use functions from the best fit
anFun_fit_dcp12_cp12 <- list("PG2_1-2-3-M0" = enr_in$anFun[[9]], #Changed ipp and adjusted to what there is in enr_in$anFun ★
                             "PG2_1-2-3-M1" = enr_in$anFun[[10]],
                             "PG2_1-2-3-M2" = enr_in$anFun[[11]],
                             "PG2_1-2-3-M3" = enr_in$anFun[[12]])

# initial parameters
#    metabolite concentrations
meta_conc_dcp12_cp12 <- c("PG2"=5.74411473678459, "PEP"=621.423909679977) #Most os Sebastian initial met concentr. 1.100000e-01 or 1.100000e-02
#    fluxes
kp_dcp12_cp12 <- c("v1"=0)
#    sd of measurements
sd_meas_dcp12_cp12 <- list(iso=sd(enr_PEP_dcp12_cp12), conc=c("PEP"=1.171822098))
#    names of parameters (concentrations & fluxes) to estimate  
te_dcp12_cp12 <- c("v1", "PEP")  # Here is somehow less parameters in the original example
#    upper & lower constraints on parameters
te_upc_dcp12_cp12 <- c("v1"=10, "PEP"=15)
te_loc_dcp12_cp12 <- c("v1"=1e-4, "PEP"=1e-1)

# run flux calculation
resF_dcp12_cp12 <- fit(net     = net,
                       times      = times,
                       kp         = kp_dcp12_cp12,
                       to_est     = te_dcp12_cp12,
                       te_upc     = te_upc_dcp12_cp12[te_dcp12_cp12],
                       te_loc     = te_loc_dcp12_cp12[te_dcp12_cp12],
                       te_upc_det = NULL,
                       te_loc_det = NULL,
                       eq_det     = NULL,
                       data_meas  = list(iso=data_exp_dcp12_cp12, conc=c("PEP"=1.171822098)),
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

# 13C-enrichments of metabolic intermediates (PEP)
enr_PEP_dcp12_cp12dCysC <- c(0, 0.1172929425, 0.00006518, 0.0141154775, 0.05690679, 0.05690679, 0.058643385)
data_exp_dcp12_cp12dCysC <- cbind(times, enr_PEP_dcp12_cp12dCysC)
colnames(data_exp_dcp12_cp12dCysC) <- c("time", "PEP_1-2-3")


# define analytical functions of label input
# here we use functions from the best fit
anFun_fit_dcp12_cp12dCysC <- list("PG2_1-2-3-M0" = enr_in$anFun[[13]], #Changed ipp and adjusted to what there is in enr_in$anFun ★
                                  "PG2_1-2-3-M1" = enr_in$anFun[[14]],
                                  "PG2_1-2-3-M2" = enr_in$anFun[[15]],
                                  "PG2_1-2-3-M3" = enr_in$anFun[[16]])

# initial parameters
#    metabolite concentrations
meta_conc_dcp12_cp12dCysC <- c("PG2"=10.611146069799, "PEP"=622.577331776105) #Most os Sebastian initial met concentr. 1.100000e-01 or 1.100000e-02
#    fluxes
kp_dcp12_cp12dCysC <- c("v1"=1) #changed to 1
#    sd of measurements
sd_meas_dcp12_cp12dCysC <- list(iso=sd(enr_PEP_dcp12_cp12dCysC), conc=c("PEP"=0.584926919))
#    names of parameters (concentrations & fluxes) to estimate  
te_dcp12_cp12dCysC <- c("v1", "PEP")  # Here is somehow less parameters in the original example
#    upper & lower constraints on parameters
te_upc_dcp12_cp12dCysC <- c("v1"=10, "PEP"=5) #changed v1 to 1.
te_loc_dcp12_cp12dCysC <- c("v1"=1e-4, "PEP"=1e-1)

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
                            data_meas  = list(iso=data_exp_dcp12_cp12dCysC, conc=c("PEP"=0.584926919)),
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
enr_PEP_dcp12_cp12dCysN <- c(0, 0, 0, 0, 0, 0.001211588, 0.001158125)
data_exp_dcp12_cp12dCysN <- cbind(times, enr_PEP_dcp12_cp12dCysN)
colnames(data_exp_dcp12_cp12dCysN) <- c("time", "PEP_1-2-3")


# define analytical functions of label input
# here we use functions from the best fit
anFun_fit_dcp12_cp12dCysN <- list("PG2_1-2-3-M0" = enr_in$anFun[[17]], #Changed ipp and adjusted to what there is in enr_in$anFun ★
                                  "PG2_1-2-3-M1" = enr_in$anFun[[18]],
                                  "PG2_1-2-3-M2" = enr_in$anFun[[19]],
                                  "PG2_1-2-3-M3" = enr_in$anFun[[20]])

# initial parameters
#    metabolite concentrations
meta_conc_dcp12_cp12dCysN <- c("PG2"=6.03444476511343, "PEP"=726.030508701045) #Most os Sebastian initial met concentr. 1.100000e-01 or 1.100000e-02
#    fluxes
kp_dcp12_cp12dCysN <- c("v1"=1)
#    sd of measurements
sd_meas_dcp12_cp12dCysN <- list(iso=sd(enr_PEP_dcp12_cp12dCysN), conc=c("PEP"=1.552505587))
#    names of parameters (concentrations & fluxes) to estimate  
te_dcp12_cp12dCysN <- c("v1", "PEP")  # Here is somehow less parameters in the original example
#    upper & lower constraints on parameters
te_upc_dcp12_cp12dCysN <- c("v1"=10, "PEP"=5)
te_loc_dcp12_cp12dCysN <- c("v1"=1e-4, "PEP"=1e-1)

# run flux calculation
resF_dcp12_cp12dCysN <- fit(net     = net,
                            times      = times,
                            kp         = kp_dcp12_cp12dCysN,
                            to_est     = te_dcp12_cp12dCysN,
                            te_upc     = te_upc_dcp12_cp12dCysN[te_dcp12_cp12dCysN],
                            te_loc     = te_loc_dcp12_cp12dCysN[te_dcp12_cp12dCysN],
                            te_upc_det = NULL,
                            te_loc_det = NULL,
                            eq_det     = NULL,
                            data_meas  = list(iso=data_exp_dcp12_cp12dCysN, conc=c("PEP"=1.552505587)),
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
enr_PEP_dcp12_cp12dCysNC <- c(0,	0,	0.002868095,	0.019125045,	0.085424458,	0.273779933,	0.265068757)
data_exp_dcp12_cp12dCysNC <- cbind(times, enr_PEP_dcp12_cp12dCysNC)
colnames(data_exp_dcp12_cp12dCysNC) <- c("time", "PEP_1-2-3")


# define analytical functions of label input
# here we use functions from the best fit
anFun_fit_dcp12_cp12dCysNC <- list("PG2_1-2-3-M0" = enr_in$anFun[[21]], #Changed ipp and adjusted to what there is in enr_in$anFun ★
                                   "PG2_1-2-3-M1" = enr_in$anFun[[22]],
                                   "PG2_1-2-3-M2" = enr_in$anFun[[23]],
                                   "PG2_1-2-3-M3" = enr_in$anFun[[24]])

# initial parameters
#    metabolite concentrations
meta_conc_dcp12_cp12dCysNC <- c("PG2"=5.50614272562167, "PEP"=728.604120288041) #Most os Sebastian initial met concentr. 1.100000e-01 or 1.100000e-02
#    fluxes
kp_dcp12_cp12dCysNC <- c("v1"=1)
#    sd of measurements
sd_meas_dcp12_cp12dCysNC <- list(iso=sd(enr_PEP_dcp12_cp12dCysNC), conc=c("PEP"=1.078894006))
#    names of parameters (concentrations & fluxes) to estimate  
te_dcp12_cp12dCysNC <- c("v1", "PEP")  # Here is somehow less parameters in the original example
#    upper & lower constraints on parameters
te_upc_dcp12_cp12dCysNC <- c("v1"=10, "PEP"=15)
te_loc_dcp12_cp12dCysNC <- c("v1"=1e-4, "PEP"=1e-1)

# run flux calculation
resF_dcp12_cp12dCysNC <- fit(net     = net,
                             times      = times,
                             kp         = kp_dcp12_cp12dCysNC,
                             to_est     = te_dcp12_cp12dCysNC,
                             te_upc     = te_upc_dcp12_cp12dCysNC[te_dcp12_cp12dCysNC],
                             te_loc     = te_loc_dcp12_cp12dCysNC[te_dcp12_cp12dCysNC],
                             te_upc_det = NULL,
                             te_loc_det = NULL,
                             eq_det     = NULL,
                             data_meas  = list(iso=data_exp_dcp12_cp12dCysNC, conc=c("PEP"=1.078894006)),
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


#load("C:/Users/coron/Documents/Thesis/isosim/models/thesis_pathway/Run-2pg-pep/Adaptation_2PG_PEP.RData")

####################################
cat("\n   ... plot results ...\n\n")
####################################

# colors for "WT", "dcp12", "dcp12::cp12", "dcp12::cp12dCysC", "dcp12::cp12dCysN and "dcp12::cp12dCysNC" datasets
col_f <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b")

# pdf(paste("Fig_2PG_PEP_WT_dcp12_dcp12_cp12_",I,".pdf", sep=""))
# par(mfrow = c(3,2))

# Small offset to avoid plotting issues with zero-only values
epsilon <- 1e-6

# Standard deviation for WT
sd_wt <- list("2PG"=c(0, 0, 0, 0, 0, 0, 0),
              "PEP"=c(0, 0, 0, 0, 0, 0, 0))
# 
# # Main plot loop
#   met <- ("PEP")
# 
#   # Add small offset to zero values for plotting
#   data_wt_offset <- data_exp_wt[,2] + epsilon
#   sd_wt_offset <- sd_wt[[met]] + epsilon
# 
#   # Plot the experimental data
#   plot(times, data_wt_offset, type="p", ylim=c(0,1), pch=21, col=col_f[1], bg=col_f[1], main=met, las=1, xlab="Time (min)", ylab="13C-Enrichment")
# 
# 
#   # Plot error bars only if sd_wt is non-zero
#   if (any(sd_wt[[met]] != 0)) {
#     segments(x0=times, y0=data_wt_offset-sd_wt_offset, x1=times, y1=data_wt_offset+sd_wt_offset)
#     segments(x0=times-1.5, y0=data_wt_offset+sd_wt_offset, x1=times+1.5, y1=data_wt_offset+sd_wt_offset)
#     segments(x0=times-1.5, y0=data_wt_offset-sd_wt_offset, x1=times+1.5, y1=data_wt_offset-sd_wt_offset)
#   }
# 
#   # Plot data points and simulation line
#   points(times, data_wt_offset, pch=21, col=col_f[1], bg=col_f[1])
#   lines(times, resF_wt$result$retres$sim + epsilon, col=col_f[1])
# 
# # Additional plot for "2PG"
# met = "2PG"
# data_2pg_offset <- data_input[,"WTK_PG2_1-2-3"] + epsilon
# sd_2pg_dcp12_offset <- sd_wt[[met]] + epsilon
# 
# # Plot data for "2PG"
# plot(times, data_2pg_offset, type="p", ylim=c(0,1), pch=21, col=col_f[1], bg=col_f[1], main=met, las=1, xlab="Time (min)", ylab="13C-Enrichment")
# 
# # Plot error bars for "2PG" if sd_wt is non-zero
# if (any(sd_wt[[met]] != 0)) {
#   segments(x0=times, y0=data_2pg_offset-sd_2pg_dcp12_offset, x1=times, y1=data_2pg_offset+sd_2pg_dcp12_offset)
#   segments(x0=times-1.5, y0=data_2pg_offset+sd_2pg_dcp12_offset, x1=times+1.5, y1=data_2pg_offset+sd_2pg_dcp12_offset)
#   segments(x0=times-1.5, y0=data_2pg_offset-sd_2pg_dcp12_offset, x1=times+1.5, y1=data_2pg_offset-sd_2pg_dcp12_offset)
# }
# 
# # Plot data points and evaluation function line
# points(times, data_2pg_offset, pch=21, col=col_f[1], bg=col_f[1])
# 
# # Evaluate and plot fitted function
# eval(parse(text=paste("fe=function(t){", anFun_fit_wt[[4]], "}", sep="")))
# lines(times, fe(times) + epsilon, col=col_f[1])
# 
# # Close PDF device
# #dev.off()
# 
# 
# #################################################################################
# 
# # dcp12
# 
# #png(paste("Fig_All_dcp12.pdf", sep=""), height=3, width=9)
# #par(mfrow = c(1,3))
# 
# # Small offset to avoid plotting issues with zero-only values
# epsilon <- 1e-6

# Standard deviation for dcp12
sd_dcp12 <- list("2PG"=c(0, 0, 0.002047375, 0.007319948, 0.034679842, 0.043441443, 0.075407513),
                 "PEP"=c(0, 0, 0.011553127, 0.012540274, 0.020043716, 0.043465652, 0.030065502)
)

# # Main plot loop for dcp12
#   met <- c("PEP")
# 
#   # Add small offset to zero values for plotting
#   data_dcp12_offset <- data_exp_dcp12[,2] + epsilon
#   sd_dcp12_offset <- sd_dcp12[[met]] + epsilon
# 
#   # Plot the experimental data with x and y labels
#   plot(times, data_dcp12_offset, type="p", ylim=c(0,1), pch=21, col=col_f[2], bg=col_f[2], main=met, las=1,
#        xlab="Time (min)", ylab="13C-Enrichment")
# 
# 
#   # Plot error bars only if the standard deviation is non-zero
#   if (any(sd_dcp12[[met]] != 0)) {
#     segments(x0=times, y0=data_dcp12_offset-sd_dcp12_offset, x1=times, y1=data_dcp12_offset+sd_dcp12_offset)
#     segments(x0=times-1.5, y0=data_dcp12_offset+sd_dcp12_offset, x1=times+1.5, y1=data_dcp12_offset+sd_dcp12_offset)
#     segments(x0=times-1.5, y0=data_dcp12_offset-sd_dcp12_offset, x1=times+1.5, y1=data_dcp12_offset-sd_dcp12_offset)
#   }
# 
#   # Plot data points and simulation line
#   points(times, data_dcp12_offset, pch=21, col=col_f[2], bg=col_f[2])
#   lines(times, resF_dcp12$result$retres$sim + epsilon, col=col_f[2])
# 
# 
# # Additional plot for "2PG"
# met = "2PG"
# data_2pg_dcp12_offset <- data_input[,"dcp12_PG2_1-2-3"] + epsilon
# sd_2pg_dcp12_offset <- sd_dcp12[[met]] + epsilon
# 
# # Plot data for "2PG"
# plot(times, data_2pg_dcp12_offset, type="p", ylim=c(0,1), pch=21, col=col_f[2], bg=col_f[2], main=met, las=1,
#      xlab="Time (min)", ylab="13C-Enrichment")
# 
# # Plot error bars for "2PG" if sd_dcp12 is non-zero
# if (any(sd_dcp12[[met]] != 0)) {
#   segments(x0=times, y0=data_2pg_dcp12_offset-sd_2pg_dcp12_offset, x1=times, y1=data_2pg_dcp12_offset+sd_2pg_dcp12_offset)
#   segments(x0=times-1.5, y0=data_2pg_dcp12_offset+sd_2pg_dcp12_offset, x1=times+1.5, y1=data_2pg_dcp12_offset+sd_2pg_dcp12_offset)
#   segments(x0=times-1.5, y0=data_2pg_dcp12_offset-sd_2pg_dcp12_offset, x1=times+1.5, y1=data_2pg_dcp12_offset-sd_2pg_dcp12_offset)
# }
# 
# # Plot data points
# points(times, data_2pg_dcp12_offset, pch=21, col=col_f[2], bg=col_f[2])
# 
# # Evaluate and plot fitted function
# eval(parse(text=paste("fe=function(t){", anFun_fit_dcp12[[4]], "}", sep="")))
# lines(times, fe(times) + epsilon, col=col_f[2])
# 
# # Close PDF device
# #dev.off()
# 
# ################################################################################################
# # dcp1::cp12
# 
# #png(paste("Fig_All_dcp12_cp12.pdf", sep=""), height=3, width=9)
# #par(mfrow = c(1,3))
# 
# # Small offset to avoid plotting issues with zero-only values
# epsilon <- 1e-6

# Standard deviation for dcp12::cp12
sd_dcp12_cp12 <- list(
  "2PG"=c(0, 0,	0, 0,	0, 0.000506475,	0.00610141),
  "PEP"=c(0, 0, 0, 0, 0, 0.00942933, 0.012240472)
)

# # Main plot loop for dcp12::cp12
#   met <- c("PEP")
# 
#   # Add small offset to zero values for plotting
#   data_dcp12_cp12_offset <- data_exp_dcp12_cp12[,2] + epsilon
#   sd_dcp12_cp12_offset <- sd_dcp12_cp12[[met]] + epsilon
# 
#   # Plot the experimental data with x and y labels
#   plot(times, data_dcp12_cp12_offset, type="p", ylim=c(0,1), pch=21, col=col_f[3], bg=col_f[3], main=met, las=1,
#        xlab="Time (min)", ylab="13C-Enrichment")
# 
# 
#   # Plot error bars only if the standard deviation is non-zero
#   if (any(sd_dcp12_cp12[[met]] != 0)) {
#     segments(x0=times, y0=data_dcp12_cp12_offset-sd_dcp12_cp12_offset, x1=times, y1=data_dcp12_cp12_offset+sd_dcp12_cp12_offset)
#     segments(x0=times-1.5, y0=data_dcp12_cp12_offset-sd_dcp12_cp12_offset, x1=times+1.5, y1=data_dcp12_cp12_offset+sd_dcp12_cp12_offset)
#     segments(x0=times-1.5, y0=data_dcp12_cp12_offset-sd_dcp12_cp12_offset, x1=times+1.5, y1=data_dcp12_cp12_offset-sd_dcp12_cp12_offset)
#   }
# 
#   # Plot data points and simulation line
#   points(times, data_dcp12_cp12_offset, pch=21, col=col_f[3], bg=col_f[3])
#   lines(times, resF_dcp12_cp12$result$retres$sim + epsilon, col=col_f[3])
# 
# 
# # Additional plot for "2PG"
# met = "2PG"
# data_2pg_dcp12_cp12_offset <- data_input[,"dcp12cp12_PG2_1-2-3"] + epsilon
# sd_2pg_dcp12_cp12_offset <- sd_dcp12_cp12[[met]] + epsilon
# 
# # Plot data for "2PG"
# plot(times, data_2pg_dcp12_cp12_offset, type="p", ylim=c(0,1), pch=21, col=col_f[3], bg=col_f[3], main=met, las=1,
#      xlab="Time (min)", ylab="13C-Enrichment")
# 
# # Plot error bars for "2PG" if sd_dcp12 is non-zero
# if (any(sd_dcp12[[met]] != 0)) {
#   segments(x0=times, y0=data_2pg_dcp12_cp12_offset-sd_2pg_dcp12_cp12_offset, x1=times, y1=data_2pg_dcp12_cp12_offset+sd_2pg_dcp12_cp12_offset)
#   segments(x0=times-1.5, y0=data_2pg_dcp12_cp12_offset+sd_2pg_dcp12_cp12_offset, x1=times+1.5, y1=data_2pg_dcp12_cp12_offset+sd_2pg_dcp12_cp12_offset)
#   segments(x0=times-1.5, y0=data_2pg_dcp12_cp12_offset-sd_2pg_dcp12_cp12_offset, x1=times+1.5, y1=data_2pg_dcp12_cp12_offset-sd_2pg_dcp12_cp12_offset)
# }
# 
# # Plot data points
# points(times, data_2pg_dcp12_cp12_offset, pch=21, col=col_f[3], bg=col_f[3])
# 
# # Evaluate and plot fitted function
# eval(parse(text=paste("fe=function(t){", anFun_fit_dcp12_cp12[[4]], "}", sep="")))
# lines(times, fe(times) + epsilon, col=col_f[3])
# 
# # Close PDF device
# dev.off()
# 
# ##############################################################################################
# # dcp12::cp12dCysC
# 
# pdf(paste("Fig_2PG_PEP_dcp12_cp12dCysC_",I,".pdf", sep=""))
# #par(mfrow = c(1,3))
# par(mfrow = c(3,2))
# 
# # Small offset to avoid plotting issues with zero-only values
# epsilon <- 1e-6

# Standard deviation for dcp12::cp12
sd_dcp12_cp12dCysC <- list(
  "2PG"=c(0, 0, 0, 0.007980192, 0.01988221, 0.011417744, 0.013267179),
  "PEP"=c(0, 0.234585885, 0.00013036, 0.006763073, 0.011860306, 0.008635114, 0.004884389)
)

# # Main plot loop for dcp12::cp12
#   met <- c("PEP")
# 
#   # Add small offset to zero values for plotting
#   data_dcp12_cp12dCysC_offset <- data_exp_dcp12_cp12dCysC[,2] + epsilon
#   sd_dcp12_cp12dCysC_offset <- sd_dcp12_cp12dCysC[[met]] + epsilon
# 
#   # Plot the experimental data with x and y labels
#   plot(times, data_dcp12_cp12dCysC_offset, type="p", ylim=c(0,1), pch=21, col=col_f[4], bg=col_f[4], main=met, las=1,
#        xlab="Time (min)", ylab="13C-Enrichment")
# 
#   # Plot error bars only if the standard deviation is non-zero
#   if (any(sd_dcp12_cp12dCysC[[met]] != 0)) {
#     segments(x0=times, y0=data_dcp12_cp12dCysC_offset-sd_dcp12_cp12dCysC_offset, x1=times, y1=data_dcp12_cp12dCysC_offset+sd_dcp12_cp12dCysC_offset)
#     segments(x0=times-1.5, y0=data_dcp12_cp12dCysC_offset-sd_dcp12_cp12dCysC_offset, x1=times+1.5, y1=data_dcp12_cp12dCysC_offset+sd_dcp12_cp12dCysC_offset)
#     segments(x0=times-1.5, y0=data_dcp12_cp12dCysC_offset-sd_dcp12_cp12dCysC_offset, x1=times+1.5, y1=data_dcp12_cp12dCysC_offset-sd_dcp12_cp12dCysC_offset)
#   }
# 
#   # Plot data points and simulation line
#   points(times, data_dcp12_cp12dCysC_offset, pch=21, col=col_f[4], bg=col_f[4])
#   lines(times, resF_dcp12_cp12dCysC$result$retres$sim + epsilon, col=col_f[4])
# 
# 
# # Additional plot for "2PG"
# met = "2PG"
# data_2pg_dcp12_cp12dCysC_offset <- data_input[,"dcp12cp12dCysC_PG2_1-2-3"] + epsilon
# sd_2pg_dcp12_cp12dCysC_offset <- sd_dcp12_cp12dCysC[[met]] + epsilon
# 
# # Plot data for "2PG"
# plot(times, data_2pg_dcp12_cp12dCysC_offset, type="p", ylim=c(0,1), pch=21, col=col_f[4], bg=col_f[4], main=met, las=1,
#      xlab="Time (min)", ylab="13C-Enrichment")
# 
# # Plot error bars for "2PG" if sd_dcp12 is non-zero
# if (any(sd_dcp12[[met]] != 0)) {
#   segments(x0=times, y0=data_2pg_dcp12_cp12dCysC_offset-sd_2pg_dcp12_cp12dCysC_offset, x1=times, y1=data_2pg_dcp12_cp12dCysC_offset+sd_2pg_dcp12_cp12dCysC_offset)
#   segments(x0=times-1.5, y0=data_2pg_dcp12_cp12dCysC_offset+sd_2pg_dcp12_cp12dCysC_offset, x1=times+1.5, y1=data_2pg_dcp12_cp12dCysC_offset+sd_2pg_dcp12_cp12dCysC_offset)
#   segments(x0=times-1.5, y0=data_2pg_dcp12_cp12dCysC_offset-sd_2pg_dcp12_cp12dCysC_offset, x1=times+1.5, y1=data_2pg_dcp12_cp12dCysC_offset-sd_2pg_dcp12_cp12dCysC_offset)
# }
# 
# # Plot data points
# points(times, data_2pg_dcp12_cp12dCysC_offset, pch=21, col=col_f[4], bg=col_f[4])
# 
# # Evaluate and plot fitted function
# eval(parse(text=paste("fe=function(t){", anFun_fit_dcp12_cp12dCysC[[4]], "}", sep="")))
# lines(times, fe(times) + epsilon, col=col_f[4])
# 
# # Close PDF device
# #dev.off()
# 
# 
# ##############################################################################################
# # dcp12::cp12dCysN
# 
# #png(paste("Fig_All_dcp12_cp12dCysN.pdf", sep=""), height=3, width=9)
# #par(mfrow = c(1,3))
# 
# # Small offset to avoid plotting issues with zero-only values
# epsilon <- 1e-6
# 
# Standard deviation for dcp12::cp12
sd_dcp12_cp12dCysN <- list(
  "2PG"=c(0, 0, 0, 0, 0, 0, 0),
  "PEP"=c(0, 0, 0, 0, 0, 0.00167140000030663, 0.00231625)
)

# # Main plot loop for dcp12::cp12
#   met <- c("PEP")
# 
#   # Add small offset to zero values for plotting
#   data_dcp12_cp12dCysN_offset <- data_exp_dcp12_cp12dCysN[,2] + epsilon
#   sd_dcp12_cp12dCysN_offset <- sd_dcp12_cp12dCysN[[met]] + epsilon
# 
#   # Plot the experimental data with x and y labels
#   plot(times, data_dcp12_cp12dCysN_offset, type="p", ylim=c(0,1), pch=21, col=col_f[5], bg=col_f[5], main=met, las=1,
#        xlab="Time (min)", ylab="13C-Enrichment")
# 
#   # Plot error bars only if the standard deviation is non-zero
#   if (any(sd_dcp12_cp12dCysN[[met]] != 0)) {
#     segments(x0=times, y0=data_dcp12_cp12dCysN_offset-sd_dcp12_cp12dCysN_offset, x1=times, y1=data_dcp12_cp12dCysN_offset+sd_dcp12_cp12dCysN_offset)
#     segments(x0=times-1.5, y0=data_dcp12_cp12dCysN_offset-sd_dcp12_cp12dCysN_offset, x1=times+1.5, y1=data_dcp12_cp12dCysN_offset+sd_dcp12_cp12dCysN_offset)
#     segments(x0=times-1.5, y0=data_dcp12_cp12dCysN_offset-sd_dcp12_cp12dCysN_offset, x1=times+1.5, y1=data_dcp12_cp12dCysN_offset-sd_dcp12_cp12dCysN_offset)
#   }
# 
#   # Plot data points and simulation line
#   points(times, data_dcp12_cp12dCysN_offset, pch=21, col=col_f[5], bg=col_f[5])
#   lines(times, resF_dcp12_cp12dCysN$result$retres$sim + epsilon, col=col_f[5])
# 
# 
# # Additional plot for "2PG"
# met = "2PG"
# data_2pg_dcp12_cp12dCysN_offset <- data_input[,"dcp12cp12dCysN_PG2_1-2-3"] + epsilon
# sd_2pg_dcp12_cp12dCysN_offset <- sd_dcp12_cp12dCysN[[met]] + epsilon
# 
# # Plot data for "2PG"
# plot(times, data_2pg_dcp12_cp12dCysN_offset, type="p", ylim=c(0,1), pch=21, col=col_f[5], bg=col_f[5], main=met, las=1,
#      xlab="Time (min)", ylab="13C-Enrichment")
# 
# # Plot error bars for "2PG" if sd_dcp12 is non-zero
# if (any(sd_dcp12[[met]] != 0)) {
#   segments(x0=times, y0=data_2pg_dcp12_cp12dCysN_offset-sd_2pg_dcp12_cp12dCysN_offset, x1=times, y1=data_2pg_dcp12_cp12dCysN_offset+sd_2pg_dcp12_cp12dCysN_offset)
#   segments(x0=times-1.5, y0=data_2pg_dcp12_cp12dCysN_offset+sd_2pg_dcp12_cp12dCysN_offset, x1=times+1.5, y1=data_2pg_dcp12_cp12dCysN_offset+sd_2pg_dcp12_cp12dCysN_offset)
#   segments(x0=times-1.5, y0=data_2pg_dcp12_cp12dCysN_offset-sd_2pg_dcp12_cp12dCysN_offset, x1=times+1.5, y1=data_2pg_dcp12_cp12dCysN_offset-sd_2pg_dcp12_cp12dCysN_offset)
# }
# 
# # Plot data points
# points(times, data_2pg_dcp12_cp12dCysN_offset, pch=21, col=col_f[5], bg=col_f[5])
# 
# # Evaluate and plot fitted function
# eval(parse(text=paste("fe=function(t){", anFun_fit_dcp12_cp12dCysN[[4]], "}", sep="")))
# lines(times, fe(times) + epsilon, col=col_f[5])
# 
# # Close PDF device
# #dev.off()
# 
# 
# 
# ##############################################################################################
# # dcp12::cp12dCysNC
# 
# #png(paste("Fig_All_dcp12_cp12dCysNC.pdf", sep=""), height=3, width=9)
# #par(mfrow = c(1,3))
# 
# # Small offset to avoid plotting issues with zero-only values
# epsilon <- 1e-6

# Standard deviation for dcp12::cp12
sd_dcp12_cp12dCysNC <- list(
  "2PG"=c(0, 0, 0, 0.019012441, 0.029334965, 0.023493797, 0.037075651),
  "PEP"=c(0, 0, 0.00573619, 0.019240409, 0.025263831, 0.004137565, 0.029436832)
)
# 
# # Main plot loop for dcp12::cp12
#   met <- c("PEP")
# 
#   # Add small offset to zero values for plotting
#   data_dcp12_cp12dCysNC_offset <- data_exp_dcp12_cp12dCysNC[,2] + epsilon
#   sd_dcp12_cp12dCysNC_offset <- sd_dcp12_cp12dCysNC[[met]] + epsilon
# 
#   # Plot the experimental data with x and y labels
#   plot(times, data_dcp12_cp12dCysNC_offset, type="p", ylim=c(0,1), pch=21, col=col_f[6], bg=col_f[6], main=met, las=1,
#        xlab="Time (min)", ylab="13C-Enrichment")
# 
#   # Plot error bars only if the standard deviation is non-zero
#   if (any(sd_dcp12_cp12dCysNC[[met]] != 0)) {
#     segments(x0=times, y0=data_dcp12_cp12dCysNC_offset-sd_dcp12_cp12dCysNC_offset, x1=times, y1=data_dcp12_cp12dCysNC_offset+sd_dcp12_cp12dCysNC_offset)
#     segments(x0=times-1.5, y0=data_dcp12_cp12dCysNC_offset-sd_dcp12_cp12dCysNC_offset, x1=times+1.5, y1=data_dcp12_cp12dCysNC_offset+sd_dcp12_cp12dCysNC_offset)
#     segments(x0=times-1.5, y0=data_dcp12_cp12dCysNC_offset-sd_dcp12_cp12dCysNC_offset, x1=times+1.5, y1=data_dcp12_cp12dCysNC_offset-sd_dcp12_cp12dCysNC_offset)
#   }
# 
#   # Plot data points and simulation line
#   points(times, data_dcp12_cp12dCysNC_offset, pch=21, col=col_f[6], bg=col_f[6])
#   lines(times, resF_dcp12_cp12dCysNC$result$retres$sim + epsilon, col=col_f[6])
# 
# 
# # Additional plot for "2PG"
# met = "2PG"
# data_2pg_dcp12_cp12dCysNC_offset <- data_input[,"dcp12cp12dCysNC_PG2_1-2-3"] + epsilon
# sd_2pg_dcp12_cp12dCysNC_offset <- sd_dcp12_cp12dCysNC[[met]] + epsilon
# 
# # Plot data for "2PG"
# plot(times, data_2pg_dcp12_cp12dCysNC_offset, type="p", ylim=c(0,1), pch=21, col=col_f[6], bg=col_f[6], main=met, las=1,
#      xlab="Time (min)", ylab="13C-Enrichment")
# 
# # Plot error bars for "2PG" if sd_dcp12 is non-zero
# if (any(sd_dcp12[[met]] != 0)) {
#   segments(x0=times, y0=data_2pg_dcp12_cp12dCysNC_offset-sd_2pg_dcp12_cp12dCysNC_offset, x1=times, y1=data_2pg_dcp12_cp12dCysNC_offset+sd_2pg_dcp12_cp12dCysNC_offset)
#   segments(x0=times-1.5, y0=data_2pg_dcp12_cp12dCysNC_offset+sd_2pg_dcp12_cp12dCysNC_offset, x1=times+1.5, y1=data_2pg_dcp12_cp12dCysNC_offset+sd_2pg_dcp12_cp12dCysNC_offset)
#   segments(x0=times-1.5, y0=data_2pg_dcp12_cp12dCysNC_offset-sd_2pg_dcp12_cp12dCysNC_offset, x1=times+1.5, y1=data_2pg_dcp12_cp12dCysNC_offset-sd_2pg_dcp12_cp12dCysNC_offset)
# }
# 
# # Plot data points
# points(times, data_2pg_dcp12_cp12dCysNC_offset, pch=21, col=col_f[6], bg=col_f[6])
# 
# # Evaluate and plot fitted function
# eval(parse(text=paste("fe=function(t){", anFun_fit_dcp12_cp12dCysNC[[4]], "}", sep="")))
# lines(times, fe(times) + epsilon, col=col_f[6])
# 
# # Close PDF device
# dev.off()


#####################################################################################################
# par(mar = c(10,5,4,2))
# pdf(paste("Fig_Estimated_fluxes_2pg_pep_",I,".pdf", sep=""))

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

# barCenters <- barplot(flux, beside=FALSE, ylim=c(0,10), ylab="2PG production (nmol/OD750/min)", 
#                       names.arg=c("WT", "dcp12", "dcp12:cp12", "dcp12:cp12dCysC", "dcp12:cp12dCysN",
#                                   "dcp12:cp12dCysNC"), las=1, col=col_f, main = "Estimated fluxes 2PG-PEP")
# segments(barCenters, flux - sd_flux, barCenters,
#          flux + sd_flux, lwd = 1.5)
# segments(barCenters-0.1, flux + sd_flux, barCenters+0.1,
#          flux + sd_flux, lwd = 1.5)
# segments(barCenters-0.1, flux - sd_flux, barCenters+0.1,
#          flux - sd_flux, lwd = 1.5)
# text(x = barCenters, y = flux + sd_flux + 0.05, labels = round(flux, 4), xpd = TRUE, cex = 0.9)
# 
# dev.off()
# 
# par(mar = c(10,5,4,2))
# pdf(paste("Fig_PEP_turnover_2pg_pep_",I,".pdf", sep=""))
# 
turnover_PEP <- c(resF_wt$result$par["PEP"]/resF_wt$result$par["v1"],
                  resF_dcp12$result$par["PEP"]/resF_dcp12$result$par["v1"],
                  resF_dcp12_cp12$result$par["PEP"]/resF_dcp12_cp12$result$par["v1"],
                  resF_dcp12_cp12dCysC$result$par["PEP"]/resF_dcp12_cp12dCysC$result$par["v1"],
                  resF_dcp12_cp12dCysN$result$par["PEP"]/resF_dcp12_cp12dCysN$result$par["v1"],
                  resF_dcp12_cp12dCysNC$result$par["PEP"]/resF_dcp12_cp12dCysNC$result$par["v1"])

sd_turnover_PEP <- c(sqrt(resF_wt$sens$summary["v1", "sd"]**2+resF_wt$sens$summary["PEP", "sd"]**2),
                     sqrt(resF_dcp12$sens$summary["v1", "sd"]**2+resF_dcp12$sens$summary["PEP", "sd"]**2),
                     sqrt(resF_dcp12_cp12$sens$summary["v1", "sd"]**2+resF_dcp12_cp12$sens$summary["PEP", "sd"]**2),
                     sqrt(resF_dcp12_cp12dCysC$sens$summary["v1", "sd"]**2+resF_dcp12_cp12dCysC$sens$summary["PEP", "sd"]**2),
                     sqrt(resF_dcp12_cp12dCysN$sens$summary["v1", "sd"]**2+resF_dcp12_cp12dCysN$sens$summary["PEP", "sd"]**2),
                     sqrt(resF_dcp12_cp12dCysNC$sens$summary["v1", "sd"]**2+resF_dcp12_cp12dCysNC$sens$summary["PEP", "sd"]**2))

# barCenters <- barplot(turnover_PEP, beside=TRUE, ylim=c(0,100), ylab="PEP turnover (min-1)",
#                       names.arg=c("WT", "dcp12", "dcp12:cp12", "dcp12:cp12dCysC", "dcp12:cp12dCysN", "dcp12:cp12dCysNC"), las=1, col=col_f)
# segments(barCenters, turnover_PEP - sd_turnover_PEP, barCenters,
#          turnover_PEP + sd_turnover_PEP, lwd = 1.5)
# segments(barCenters-0.1, turnover_PEP + sd_turnover_PEP, barCenters+0.1,
#          turnover_PEP + sd_turnover_PEP, lwd = 1.5)
# segments(barCenters-0.1, turnover_PEP - sd_turnover_PEP, barCenters+0.1,
#          turnover_PEP - sd_turnover_PEP, lwd = 1.5)
# 
# # Add the turnover rate values on top of the bars
# text(x = barCenters, y = turnover_PEP + sd_turnover_PEP + 5, labels = round(turnover_PEP, 3), xpd = TRUE, cex = 0.9)
# dev.off()

####################################
cat("\nDone.\n")
####################################



