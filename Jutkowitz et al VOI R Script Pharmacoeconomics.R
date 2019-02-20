######## The Curve of Optimal Sample Size (COSS) ###########

# Developed by
# Eric Jutkowitz, PhD (1)	
# Fernando Alarid-Escudero, PhD (2) 
# Karen M. Kuntz, ScD (3)
# Hawre J. Jalal, MD, PhD (4) 

# Institutions: 		
# 1 Brown University School of Public Health, Providence RI, USA
# 2 Drug Policy Program, Center for Research and Teaching in Economics (CIDE) - CONACyT, Aguascalientes, Mexico
# 3 University of Minnesota School of Public Health, Minneapolis, MN, USA
# 4 University of Pittsburgh Graduate School of Public Health, Pittsburgh, PA, USA

#### The most up to date code could be found here: 
#### https://github.com/feralaes/COSS

#------------------------------------#
#### Load packages and functions  ####
#-----------------------------------#


library(devtools)
# install_github("DARTH-git/dampack") # Install package from github
library(dampack)  # To produce CEAC, CEAF, EVPI and other plots
library(readxl)   
library(scales)
library(reshape2) 
library(dplyr)
library(gridExtra)
library(grid)
library(ggplot2)
# setwd("...")
source("Jutkowitz et al VOI R Functions Pharmacoeconomics.R")
txtsize <- 16

#------------------------#
#### Load PSA DAtaset ####
#------------------------#

m.psa <- read.csv("Gout_PSA_Data.csv", header = TRUE)
# Number of simulations
n.sim <- nrow(m.psa)
# Strategies
n.strategies     <- 3
names.strategies <- c("Allopurinol-febuxostat sequential therapy dose-escalation",
                      "Febuxostat-allopurinol sequential therapy dose-escalation",
                      "Allopurinol only dose-escalation")
### Outcomes
Outcomes <- m.psa[, 2:(n.strategies*2+1)]
# Create matrix of Effectiveness
m.e <- m.psa[, c(3, 5, 7)]
# Create matrix of Costs
m.c <- m.psa[, c(2, 4, 6)]

### Parameters
# Matrix of parameters
m.parms <- m.psa[, (n.strategies*2+2):ncol(m.psa)]

### Make PSA object for dampack
psa_coss <- make_psa_obj(cost = m.c, effectiveness = m.e, 
                         parameters = m.parms, 
                         strategies = names.strategies)

#-----------------------------#
#### Define VOI Parameters ####
#-----------------------------#

### 1. WTP Ranges 
v.wtp <- seq(5000, 150000, by = 5000)

### 2. Discount rate
disc <- c(0.03)

### 3. Decision lifetime
LT   <- 5 
time <- seq(0, LT)

### 4. Prevalent cases
prev.factor = 1 # For sens analysis
prev     <- (8.3*.30) * prev.factor   #1.6736 # In millions

### 5. Incident cases
incid.factor = 1.0  # for sens analysis
incid    <- (.30*.199) * incid.factor  #0.2*29.376e-3 # In millions

### Research design

  # Design 1 - Randomized controlled trial comparing allopurinol dose-escalation against placebo.

  # Design 2 - Observational study evaluating health utilities of gout patients.

### 6. Fixed cost of research
cost.factor <- 1
#fixed cost Design 1
fixed.cost.d1= cost.factor * 8.740032

#fixed cost Design 2
fixed.cost.d2 <- cost.factor * 0.050000

### 7. Variable cost of research
#variable cost Design 1
cost.factor.variable <- 1
cost.per.patient.d1 = cost.factor.variable * 0.008438 

#fixed cost Design 2
cost.per.patient.d2 <- cost.factor.variable *  0.000500

#-----------------------------#
### eFigure 1: CEAC & CEAF  ###
#-----------------------------#
### Generate CEAC and CEAF
gg.ceaf <- ceac(wtp = v.wtp, psa = psa_coss) # from dampack
### Plot CEAC and CEAF
gg.ceaf <- plot(gg.ceaf, n_x_ticks = 12, txtsize = 16) +
  guides(color=guide_legend(ncol=1),
         shape=guide_legend(ncol=1)) +
  theme(legend.position = "bottom")
gg.ceaf
ggsave("figs/eFig1_CEAF.jpg", gg.ceaf, width = 15, height = 11)

#---------------------------#
#### Population Analysis ####
#---------------------------#

### Total population afectd by technology calculated with `TotPop` function
tot.pop <- TotPop(time,    # Function
                  prev, 
                  incid, 
                  disc) 

#---------------------------------#
### eFigure 2: Population EVPI  ###
#---------------------------------#
df.evpi.pop <- calc_evpi(wtp = v.wtp, psa = psa_coss, tot.pop) # from dampack

my_grob <- grobTree(textGrob(paste("Discount factor: ", disc,"\n", "Technology Life Time: ",LT," years",sep=""), 
                            x = 0.05,  y = 0.9, hjust = 0,
                            gp = gpar(col = "black", fontsize = (txtsize-2), fontface = "italic")))

gg.evpi.pop <- plot(df.evpi.pop, txtsize = 16) +
  annotation_custom(my_grob)
gg.evpi.pop
ggsave("figs/eFig2_popEVPI.jpg", gg.evpi.pop, width = 15, height = 11)

#-------------------------------------------------------------#
### eFigure 3: Population EVPPI of different study designs  ###
#-------------------------------------------------------------#

### Design 1 - Randomized controlled trial comparing allopurinol dose-escalation against placebo.
## Select parameters of Allopurinol dose-escalation for EVPPI 

selParams.allo <- c(which(colnames(m.parms) == "p_allo_success_de"), 
                    which(colnames(m.parms) == "p_hyper_D"),
                    which(colnames(m.parms) == "p_die_hyper_D"),
                    which(colnames(m.parms) == "t_max_fail"),
                    which(colnames(m.parms) == "p_allo_ae"),
                    which(colnames(m.parms) == "p_stop_low"),
                    which(colnames(m.parms) == "p_stop_med"),
                    which(colnames(m.parms) == "p_stop_hi"),
                    which(colnames(m.parms) == "p_flare_uncon"),
                    which(colnames(m.parms) == "p_flare"),
                    which(colnames(m.parms) == "p_flare_uncon_off"),
                    which(colnames(m.parms) == "p_sto_ae_allo")) 
evppi.allo <- EVPPI(wtp = v.wtp,  # Function 
                    names.strategies, Outcomes, 
                    Parameters = m.parms, 
                    selParams = selParams.allo, 
                    parm = "Study Design 1: Randomized trial allopurinol dose-escalation vs placebo")
 evppi.allo.plot <- plotEVPPI(evppi.allo) # Function
 evppi.allo.plot

### Design 2 - Observational study evaluating health utilities of gout patients.
## Select parameters of health utility for EVPPI 
selParams.u <- c(which(colnames(m.parms)=="u_disutility_delta"), 
                 which(colnames(m.parms)=="u_disutility_off"),
                 which(colnames(m.parms)=="u_flare"),
                 which(colnames(m.parms)=="u_flare_ae"),
                 which(colnames(m.parms)=="u_hyper")) 
evppi.u <- EVPPI(wtp = v.wtp,  # Function 
                 names.strategies, Outcomes, 
                 Parameters = m.parms, 
                 selParams = selParams.u, 
                 parm = "Study Design 2: Observational study evaluating health utilities")
 evppi.u.plot <- plotEVPPI(evppi.u) # Function
 evppi.u.plot

## Plot EVPPI for both study designs on same graph
df.evppi <- rbind(evppi.u,
                  evppi.allo)
df.evppi.pop <- popEVPPI(df.evppi,    # Function
                          tot.pop)
gg.evppi.pop <- popEVPPIplot(df.evppi.pop,
                             disc,
                             LT, txtsize = 16) +
  scale_x_continuous(breaks = number_ticks(12)) +
  guides(shape=guide_legend(title="Data collection study design"),
         color=guide_legend(title="Data collection study design")) + 
  ggtitle("") +
  theme_bw(base_size = 16) +
  theme(legend.position = "bottom")
gg.evppi.pop
ggsave("figs/eFig3_popEVPPI.jpg", gg.evppi.pop, width = 15, height = 11)

#----------------------#
#### Figure 1: COSS ####
#----------------------#

### EVIS estimated using the method proposed by: 
###  Jalal H, Alarid-Escudero F. 
###  A Gaussian Approximation Approach for Value of Information Analysis. 
###  Med. Decis. Making. 2018;38:174-88.
### requires effective sample size used to inform the base case value 

### Effective Sample Sizes (ESS) 
n0 <- c(263, #u_disutility_delta
      2,   #t_max_fail
      64,  #p_allo_success_de
      212, #u_disutility_off
      1835,#p_hyper_D
      101, #p_die_hyper_D
      253, #p_allo_ae
      735, #p_flare
      735, #p_flare_uncon
      735, #p_flare_uncon_off
      2,   #p_stop_low
      2,   #p_stop_med
      2,   #p_stop_hi
      254, #p_stop_ae_allo
      386, #u_flare
      10, #u_flare_ae, 
      898) #u_hyper, 


### COSS for Allopurinol DE
n0.allo <- n0[selParams.allo]
## Vector with sample sizes to evaluate EVSI & COSS
n.allo <- seq(0, 5000, by = 5)  # chaning by to higher values will decrease computation time but result in a loss of percision for n*
## Cost of research
cost.res.allo <- CostRes(fixed.cost = fixed.cost.d1,
                         samp.size = n.allo,  # vector  
                         cost.per.patient = cost.per.patient.d1, 
                         delta.C = 0.001298,  # obtained from source analysis see Jutkowitz et al. Cost-Effectivness of Gout Therapies. Annals of Internal Medicine (2014)
                         delta.E = -0.867, # the cost function calculates the absolute value of the incremental net benefit. See Wilson et al. A Practical Guide to Value of Information Analysis. PharmacoEconomics (2015)
                         wtp = v.wtp/1000000, # convert wtp to units of millions 
                         clin.trial = TRUE) # 2 default arms in clincal trial default; to change default to clincal trial include argument, n.arms = 2

system.time( # This might take a few minutes
  coss.allo <- COSS(wtp = v.wtp, 
                    Strategies = names.strategies, 
                    Outcomes = Outcomes, 
                    Parameters = m.parms, 
                    tot.pop = tot.pop, 
                    selParams = selParams.allo,
                    n0 = n0.allo, 
                    samp.size = n.allo,
                    cost.res = cost.res.allo,
                    parm = "Study Design 1: Randomized trial allopurinol dose-escalation vs placebo")
)

system.time( # This might take a few minutes
  enbs.allo <- ENBS(wtp = v.wtp, 
                    Strategies = names.strategies, 
                    Outcomes = Outcomes, 
                    Parameters = m.parms, 
                    tot.pop = tot.pop, 
                    selParams = selParams.allo,
                    n0 = n0.allo, 
                    samp.size = n.allo,
                    cost.res = cost.res.allo,
                    parm = "Study Design 1: Randomized trial allopurinol dose-escalation vs placebo")
)



### COSS for Utilities
## Effective Sample size
n0.u <- n0[selParams.u]
## Vector with sample sizes to evaluate EVSI & COSS
n.u <- seq(0, 12000, by = 5) # chaning by to higher values will decrease computation time but result in a loss of percision for n*
## Cost of research
cost.res.u <- CostRes(fixed.cost.d2,
                      samp.size = n.u, 
                      cost.per.patient.d2,
                      delta.C = 0, # not an rct so no loss from putting patients in an inferior arm
                      delta.E = 0, # not an rct so no loss from putting patients in an inferior arm
                      wtp = v.wtp,
                      clin.trial = FALSE) # In Millions 1033.0522e-6
system.time( # This might take a few minutes
  coss.u <- COSS(wtp = v.wtp, 
                 Strategies = names.strategies, 
                 Outcomes = Outcomes, 
                 Parameters = m.parms, 
                 tot.pop = tot.pop, 
                 selParams = selParams.u,
                 n0 = n0.u, 
                 samp.size = n.u,
                 cost.res = cost.res.u,
                 parm = "Study Design 2: Observational study evaluating health utilities")
)

system.time( # This might take a few minutes
  enbs.u <- ENBS(wtp = v.wtp, 
                 Strategies = names.strategies, 
                 Outcomes = Outcomes, 
                 Parameters = m.parms, 
                 tot.pop = tot.pop, 
                 selParams = selParams.u,
                 n0 = n0.u, 
                 samp.size = n.u,
                 cost.res = cost.res.u,
                 parm = "Study Design 2: Observational study evaluating health utilities")
)



##Plot COSS

#Combine OSS of each study desing in dataframe for plotting
coss.aggregated <- rbind(coss.u,
                         coss.allo)

my_grob.coss <- grobTree(textGrob(paste("Discount factor: ", disc,"\n", "Technology Life Time: ",LT," years",sep=""), 
                             x = 0.1,  y = 0.9, hjust = 0,
                             gp = gpar(col = "black", fontsize = (txtsize-2), fontface = "italic")))
gg.coss <- ggplot(coss.aggregated, aes(x = WTP/1000, y = OSS, 
                                                     colour = Parameter,
                                                     shape = Parameter)) +
  annotation_custom(my_grob.coss) +
  geom_point() +
  geom_line() +
  scale_color_hue("Data collection study design", l = 50) +
  scale_shape_discrete("Data collection study design") +
  # ggtitle("Curve of Optimal Sample Size by Study Design") + 
  scale_x_continuous(breaks = number_ticks(12)) +
  scale_y_continuous(breaks = number_ticks(6), labels = comma) +
  xlab("Willingness to pay (Thousand $/QALY)") +
  ylab("Optimal Sample Size (n*)") +
  theme_bw(base_size = 16) +
  theme(legend.position = "bottom") # c(0.5, 0.5)
gg.coss
ggsave("figs/Fig1_COSS.tiff", gg.coss, width = 15, height = 11)


#---------------------------------------------#
#### Figure 2: COSS SA on Cost of Reserach ####
#---------------------------------------------#

### % +/- change in the VOI parameters
sa <- 0.5

### Cost of research
## Vary cost of research per patient
cost.res.allo.low <- CostRes(fixed.cost.d1,
                             samp.size = n.allo,  # vector 
                             cost.per.patient.d1 * (1-sa), 
                             delta.C = 0.001298,
                             delta.E = -0.867,
                             wtp = v.wtp/1000000,
                             clin.trial = TRUE) 
cost.res.allo.high <- CostRes(fixed.cost.d1,
                              samp.size = n.allo,  # vector 
                              cost.per.patient.d1 * (1+sa), 
                              delta.C = 0.001298,
                              delta.E = -0.867,
                              wtp = v.wtp/1000000, 
                              clin.trial = TRUE) 


system.time( # This might take a few minutes
  coss.allo.low <- COSS(wtp = v.wtp, 
                        Strategies = names.strategies, 
                        Outcomes = Outcomes, 
                        Parameters = m.parms, 
                        tot.pop = tot.pop, 
                        selParams = selParams.allo,
                        n0 = n0.allo, 
                        samp.size = n.allo,
                        cost.res = cost.res.allo.low,
                        parm = "Study Design 1: Randomized trial allopurinol dose-escalation vs placebo")
)

system.time( # This might take a few minutes
  coss.allo.high <- COSS(wtp = v.wtp, 
                         Strategies = names.strategies, 
                         Outcomes = Outcomes, 
                         Parameters = m.parms, 
                         tot.pop = tot.pop, 
                         selParams = selParams.allo,
                         n0 = n0.allo, 
                         samp.size = n.allo,
                         cost.res = cost.res.allo.high,
                         parm = "Study Design 1: Randomized trial allopurinol dose-escalation vs placebo")
)

coss.allo.sa <- rbind(coss.allo.low,
                      coss.allo,
                      coss.allo.high)
coss.allo.sa$Cost <- rep(c("Low", "Base case", "High"), each = nrow(coss.allo))
coss.allo.sa$Cost <- ordered(coss.allo.sa$Cost,
                             levels = c("Low", "Base case", "High"),
                             labels = c("Low (-50%)", "Base case ($8,500)", "High (+50%)"))
my_grob.coss <- grobTree(textGrob(paste("Discount factor: ", disc,"\n", "Technology Life Time: ",LT," years",sep=""), 
                                  x = 0.01,  y = 0.9, hjust = 0,
                                  gp = gpar(col = "black", fontsize = (txtsize-2), fontface = "italic")))
gg.coss.allo.sa <- ggplot(coss.allo.sa, aes(x = WTP/1000, y = OSS,
                                            shape = Cost,
                                            linetype = Cost)) +
  annotation_custom(my_grob.coss) +
  geom_point(colour = "#CB4D42", size = 2.5) +
  geom_line(colour = "#CB4D42", size = 1.1) +
  scale_linetype_manual("Cost per patient", values = c(3, 1, 2)) +
  scale_shape_manual("Cost per patient", values = c(26, 16, 26)) +
  scale_color_hue(l = 50) +
  scale_x_continuous(breaks = number_ticks(10))+
  scale_y_continuous(breaks = number_ticks(25), labels = comma) +
  xlab("Willingness to pay (Thousand $/QALY)") +
  ylab("Optimal Sample Size (n*)") +
  theme_bw(base_size = 16) +
  theme(legend.position = c(0.1, 0.5))
gg.coss.allo.sa
ggsave("figs/Fig2_COSS_SA.jpg", gg.coss.allo.sa, width = 15, height = 11)

#------------------------------#
#### eFigure 4: COSS + ENBS ####
#------------------------------#

coss.aggregated["Value"] <- "OSS"

enbs.u<-enbs.u[,-c(3,4)] 
colnames(enbs.u)[3] <- "OSS"
enbs.u["Value"] <- "ENBS"

enbs.allo<-enbs.allo[,-c(3,4)] 
colnames(enbs.allo)[3] <- "OSS"
enbs.allo["Value"] <- "ENBS"

coss.aggregated_v2 <- rbind(coss.u,
                       coss.allo)

coss.aggregated_v2["Value"] <- "OSS"
coss.aggregated_v2 <-rbind(coss.aggregated_v2, enbs.u, enbs.allo)

my_grob.coss <- grobTree(textGrob(paste("Discount factor: ", disc,"\n", "Technology Life Time: ",LT," years",sep=""), 
                                  x = 0.2,  y = 0.9, hjust = 0,
                                  gp = gpar(col = "black", fontsize = (txtsize-2), fontface = "italic")))
gg.coss_enbs <- ggplot(coss.aggregated_v2, aes(x = WTP/1000, y = OSS, 
                                       colour = Parameter,
                                       shape = Parameter)) +
  #annotation_custom(my_grob.coss) +
  facet_wrap(facets = "Value", scales = "free_y", ncol =1, 
             labeller = as_labeller(c(OSS = "Optimal Sample Size (n*)", ENBS = "Expected Net Benefit of Sampling ($ Millions)") ),
             strip.position = "left") +
  geom_point() +
  geom_line() +
  scale_color_hue("Data collection study design", l = 50) +
  scale_shape_discrete("Data collection study design") +
  # ggtitle("Curve of Optimal Sample Size by Study Design") + 
  scale_x_continuous(breaks = number_ticks(12)) +
  scale_y_continuous(breaks = number_ticks(10), labels = comma) +
  xlab("Willingness to pay (Thousand $/QALY)") +
  ylab(NULL) +
  theme_bw(base_size = 16) +
  theme(legend.position = "bottom", strip.background = element_blank(), strip.placement = "outside") # c(0.5, 0.5)
gg.coss_enbs
ggsave("figs/eFig4_COSS_ENBS.jpg", gg.coss_enbs, width = 15, height = 11)

