
#### Population estimate ####
TotPop <- function(time, prev, incid, disc = 0){
  # Computes total population afected by technology
  #
  # Args:
  #   time:  vector with time points defining technology lifetime
  #   prev:  present prevalence
  #   incid: incidence
  #   disc:  discount factor; deafult = 0.
  #
  # Returns:
  #   tot.pop: total population afected by technology over technology lifetime
  #  
  # Technology Life Time, the last entry of vector `time`
  LT            <- time[length(time)]
  # Vector with population afected by the technolgy at each time point
  pop.time      <- c(prev, rep(incid, (length(time)-1))) 
  # Vector with present value of population afected by the technolgy at each time point
  disc.pop.time <- pop.time/(1+disc)^time
  # Total population afected by the technology
  tot.pop <-sum(disc.pop.time)
}

#### Cost of Research ####
CostRes <- function(fixed.cost = 0, 
                    samp.size, 
                    cost.per.patient, 
                    delta.C,
                    delta.E,
                    wtp,
                    clin.trial = TRUE, n.arms = 2){
  # Computes the cost of collecting information (i.e., through a research study)
  #
  # Args:
  #   fixed.cost:       fixed cost of collecting information
  #                     (e.g., fixed cost of a clinical trial); default = 0
  #   samp.size:               vector with sample sizes
  #   cost.per.patient: cost per patient in research study
  #   delta.C:          Incremental Costs
  #   delta.E:          Incremental Effectiveness
  #   clin.trial:       indicator whether calculation is for a clinical trial;
  #                     default = TRUE
  #   n.arms:           Number of arms in research study design; default = 2
  #
  # Returns:
  #   cost.res: vector with the total cost of collecting information for each simple size
  #
  if (clin.trial){
    v.INMB <- abs((delta.E*wtp) - delta.C)
    Cost.Res <- fixed.cost + (n.arms*samp.size*cost.per.patient) + outer(samp.size, v.INMB)
  } else { # E.g., cohort study
    Cost.Res <- outer(fixed.cost + samp.size*cost.per.patient, rep(1, length(wtp)))
  }
  colnames(Cost.Res) <- wtp
  return(Cost.Res)
}

#### EVPPI ####
EVPPI <- function(wtp, Strategies, Outcomes, 
                  Parameters, selParams = NULL, 
                  parm = "Subset", ...){
  # Computes EVPPI for a set of parameters for different values of wtp
  #
  # Args:
  #   wtp:        Vector of willingness to pay values
  #   Strategies: Vector with strategies' names
  #   Outcomes:   Matrix with outcomes ordered in a way that for each strategy the cost must 
  #               appear first then the effectiveness
  #   Parameters: Matrix with parameters' simulations
  #   selParams:  Vector of subset of parameters to aggregate into one evppi measure.
  #               Default = NULL will perform EVPPI for each parameter separately
  #
  # Returns:
  #   evppi: Data frame with population EVPPI for different values of wtp
  # Load dependencies
  require(matrixStats)
  # Generate local variables
  nStrategy <- length(Strategies) # Number of strategies
  X         <- Parameters # Matrix with input parameters
  nParam    <- ncol(X) # Number of inputs
  n.sim      <- nrow(X) # Number of simulations
  nameParam <- colnames(X) # Name of input parameters
  
  costInd <- seq(1, 2*nStrategy, by = 2) #vector to index costs
  effInd  <- seq(2, 2*nStrategy, by = 2) #vector to index effectiveness
  if(is.null(selParams)){
    evppi <- matrix(0, nParam, length(wtp)) #vector to store EVPPI for each value of wtp
    print("Computing EVPPI for all parameters")
  } else {
    evppi <- data.frame(Parameter = parm,
                        WTP = wtp,
                        EVPPI = matrix(0, length(wtp), 1))
    print(paste("Computing EVPPI for", parm, "which consists of these parameters:", 
                nameParam[selParams]))
  }
  X0      <- matrix(1, n.sim, 1) # Intercept for LRM
  
  for(l in 1:length(wtp)){
    NMB <-  wtp[l]*Outcomes[, effInd]-Outcomes[, costInd] # Effectiveness minus Costs, with vector indexing
    ## Determine the optimal strategy with current info
    dStar <- which.max(colMeans(NMB)) #Compute E[NMB] for all strategies and determine the optimal
    ## Calculate the opportunity loss from choosing dStar
    Loss <- matrix(0, n.sim, nStrategy)
    for (s in 1:nStrategy){
      Loss[, s] <- NMB[, s] - NMB[, dStar]
    }
    if(is.null(selParams)){
      for (p in 1:nParam){
        Xp       <- cbind(X0, X[, p])
        Lhatp    <- Xp %*% (solve(t(Xp) %*% Xp) %*% t(Xp) %*% Loss) # Predict Loss using LRM
        evppi[p, l] <- mean(rowMaxs(Lhatp))
      }
    } else {
      Xp   <- data.matrix(cbind(X0, X[, selParams])) # Selected vectors
      Lhatp <- Xp %*% (solve(t(Xp) %*% Xp) %*% t(Xp) %*% Loss) # Predict Loss using LRM
      evppi[l, 3] <- mean(rowMaxs(Lhatp))
    }
  }
  if(is.null(selParams)){  
    evppi.df <- data.frame(Parameter = nameParam, evppi) # Create a dataframe with Parameters names as id var
    #colnames(evppi.df)[-1] <- paste(dollar(wtp), "/QALY", sep="") # Use it if we want the $/QALY in the wtp label
    colnames(evppi.df)[-1] <- wtp/1000
    # Calculate total EVPPI across wtp to use for sorting
    evppi.df$totEVPPI <- rowSums(evppi.df[,-1]) 
    # Order EVPPI in decreasing order
    evppi.df$Parameter <- factor(as.character(evppi.df$Parameter), 
                                 levels = evppi.df$Parameter[order(evppi.df$totEVPPI, 
                                                                   decreasing = TRUE)])
    # melt EVPPI.df but the last column with total EVPPI
    evppi.df.lng <- melt(evppi.df[, -ncol(evppi.df)] , id.vars = "Parameter") 
    colnames(evppi.df.lng)[-1] <- c("WTP", "EVPPI")
    return(evppi.df.lng)
  } else {
    evppi$WTP <- evppi$WTP/1000
    return(evppi)  
  }
}

plotEVPPI<-function(evppi,
                    n.ticks.x = 5,
                    n.ticks.y = 6,
                    title = "Expected Value of Partial Perfect Information", 
                    txtsize = 12){
  # Plots an EVPPI data frame with an EVPPI for each paramater or 
  # subset of parameters at each willingness to pay. The function 
  # differenatiates between a data frame with only one parameter and more than
  # one parameter
  # Args:
  #   evppi:      data frame with an EVPPI of each paramater at each 
  #               willingness to pay
  # Returns:
  #   ggplot object
  if(nlevels(evppi$Parameter) > 1){
    ggplot(data = evppi, aes(x = WTP, y = EVPPI, 
                             group = Parameter, 
                             colour = Parameter,
                             #linetype = factor(Parameter),
                             shape = Parameter)) +
      geom_point() +
      geom_line() +
      ggtitle(title) + 
      scale_x_discrete(breaks = as.factor(seq(25, 400, by = 25)), labels = comma)+
      scale_y_continuous(breaks = number_ticks(n.ticks.y), labels = comma) +
      scale_colour_manual(values = 1:nlevels(evppi$Parameter)) +
      scale_shape_manual(values = 1:nlevels(evppi$Parameter)) +
      scale_linetype_manual(values = 1:nlevels(evppi$Parameter)) +
      xlab("Willingness to pay (Thousand $/QALY)") +
      ylab("EVPPI ($)") +
      theme_bw() +
      theme(legend.title = element_text(size = txtsize),
            legend.key = element_rect(colour = "black"),
            legend.text = element_text(size = txtsize),
            title = element_text(face="bold", size=15),
            axis.title.x = element_text(face="bold", size=txtsize),
            axis.title.y = element_text(face="bold", size=txtsize),
            axis.text.y = element_text(size=txtsize),
            axis.text.x = element_text(size=txtsize))
  } else{
    ggplot(data = evppi, aes(x = WTP, y = EVPPI, 
                             group = Parameter, 
                             colour = Parameter,
                             #linetype = factor(Parameter),
                             shape = Parameter)) +
      geom_point() +
      geom_line() +
      ggtitle(title) + 
      scale_x_discrete(breaks = as.factor(seq(25, 400, by = 25)), labels = comma)+
      scale_y_continuous(breaks = number_ticks(n.ticks.y), labels = comma) +
      xlab("Willingness to pay (Thousand $/QALY)") +
      ylab("EVPPI ($)") +
      theme_bw() +
      theme(legend.title = element_text(size = txtsize),
            legend.key = element_rect(colour = "black"),
            legend.text = element_text(size = txtsize),
            title = element_text(face="bold", size=15),
            axis.title.x = element_text(face="bold", size=txtsize),
            axis.title.y = element_text(face="bold", size=txtsize),
            axis.text.y = element_text(size=txtsize),
            axis.text.x = element_text(size=txtsize))
  }
}


popEVPPI <- function(evppi, tot.pop){
  # Computes population EVPPI
  #
  # Args:
  #   evppi:   Data frame with EVPPI for different values of lambdas (from `EVPPI` function)
  #   tot.pop: total population afected by technology over technology lifetime
  #
  # Returns:
  #   popEVPPI: Data frame with population EVPPI for different values of lambdas
  #
  pop.evppi <- evppi
  pop.evppi$popEVPPI <- pop.evppi$EVPPI*tot.pop
  #pop.evppi <- pop.evppi[,-3]
  
  return(pop.evppi)
}

popEVPPIplot <- function(pop.evppi, disc, LT, 
                       title = "Expected Value of Partial Perfect Information", txtsize = 12){
  # Function to plot a population EVPPI data frame
  # EVPPI for each willingness to pay (lambda)
  my_grob = grobTree(textGrob(paste("Discount factor: ", disc,"\n", "Technology Life Time: ",LT," years",sep=""), 
                              x = 0.05,  y = 0.9, hjust = 0,
                              gp = gpar(col = "black", fontsize = (txtsize-2), fontface = "italic")))
  if(nlevels(pop.evppi$Parameter) > 1){
  ggplot(data = pop.evppi, aes(x = WTP, y = popEVPPI, 
                           group = Parameter,colour = Parameter,
                           #linetype = factor(Parameter),
                           shape = Parameter)) +
    geom_point() +
    geom_line() +
    ggtitle(title) + 
    annotation_custom(my_grob) +
    scale_x_discrete(breaks = as.factor(seq(25, 400, by = 25)), labels = comma)+
    scale_y_continuous(breaks = number_ticks(6), labels = comma)+
    #scale_colour_discrete("Parameter", l=50) +
    scale_colour_manual(values = 1:nlevels(pop.evppi$Parameter)) +
    scale_shape_manual(values = 1:nlevels(pop.evppi$Parameter)) +
    scale_linetype_manual(values = 1:nlevels(pop.evppi$Parameter)) +
    xlab("Willingness to pay (Thousand $/QALY)") +
    ylab("Population EVPPI (Million $)") +
    theme_bw() +
    theme(legend.title = element_text(size = txtsize),
          legend.key = element_rect(colour = "black"),
          legend.text = element_text(size = txtsize),
          title = element_text(face="bold", size=15),
          axis.title.x = element_text(face="bold", size=txtsize),
          axis.title.y = element_text(face="bold", size=txtsize),
          axis.text.y = element_text(size=txtsize),
          axis.text.x = element_text(size=txtsize))
  } else {
    ggplot(data = pop.evppi, aes(x = WTP, y = popEVPPI, 
                             group = Parameter,colour = Parameter,
                             #linetype = factor(Parameter),
                             shape = Parameter)) +
      geom_point() +
      geom_line() +
      ggtitle(title) + 
      annotation_custom(my_grob) +
      scale_x_discrete(breaks = as.factor(seq(25, 400, by = 25)), labels = comma)+
      scale_y_continuous(breaks = number_ticks(6), labels = comma)+
      xlab("Willingness to pay (thousand $/QALY)") +
      ylab("Population EVPPI (Million $)") +
      theme_bw() +
      theme(legend.title = element_text(size = txtsize),
            legend.key = element_rect(colour = "black"),
            legend.text = element_text(size = txtsize),
            title = element_text(face="bold", size=15),
            axis.title.x = element_text(face="bold", size=txtsize),
            axis.title.y = element_text(face="bold", size=txtsize),
            axis.text.y = element_text(size=txtsize),
            axis.text.x = element_text(size=txtsize))    
  }
}

#### EVSI ####
EVPSI <- function(wtp, Strategies, Outcomes, Parameters, 
                  selParams, n0, samp.size, ...){
  # Computes EVPSI for a set of parameters an specific research studies
  # sample sizes
  #
  # Args:
  #   wtp:        Vector of willingness to pay values
  #   Strategies: Vector with strategies' names
  #   Outcomes:   Matrix with outcomes ordered in a way that for each strategy the cost must 
  #               appear first then the effectiveness
  #   Parameters: Matrix with parameters' simulations
  #   selParams:  Vector of subset of parameters to aggregate into one evpsi measure 
  #               for each WTP at each samp.size
  #   n0:         Vector with effective sample size of parameters
  #   samp.size:  Research study sample sizes
  #
  # Returns:
  #   evpsi: Data frame with population EVPSI for different values of WTP
  #
  nStrategy <- length(Strategies) # Number of strategies
  X         <- Parameters # Matrix with input parameters
  nParam    <- ncol(X) # Number of inputs
  nSim      <- nrow(X) # Number of simulations
  nameParam <- colnames(X) # Name of input parameters
  nSample   <- length(samp.size)  # Number of samples
  nWTP      <- length(wtp)  # Number of WTP
  
  costInd <- seq(1, 2*nStrategy, by = 2) #vector to index costs
  effInd  <- seq(2, 2*nStrategy, by = 2) #vector to index effectiveness
  evpsi   <- array(0, dim = c(nSample, nWTP),
                   dimnames = list(samp.size, wtp))
  X0 <- matrix(1, nSim, 1) # Intercept for LRM
  Xp <- data.matrix(cbind(X0, X[, selParams])) # Selected vectors
  
  for(l in 1:nWTP){
    NMB <-  wtp[l]*Outcomes[, effInd]-Outcomes[, costInd] # Effectiveness minus Costs, with vector indexing
    ## Determine the optimal strategy with current info
    dStar <- which.max(colMeans(NMB)) #Compute E[NMB] for all strategies and determine the optimal
    ## Calculate the opportunity loss from choosing dStar
    Loss <- matrix(0, nSim, nStrategy)
    for (s in 1:nStrategy){
      Loss[, s] <- NMB[, s] - NMB[, dStar]
    }
    # Mean of selected vectors
    Xpmean <- colMeans(Xp[, -1])  
    Bp <- solve(t(Xp) %*% Xp) %*% t(Xp) %*% Loss  # Regression coefficient matrix
    for (nSamp in 1:nSample){
      v      <- samp.size[nSamp] / (samp.size[nSamp] + n0)
      sqrtv  <- sqrt(v)
      # Shrink each vector of parameters by the corresponding element in v
      #Xp.1 <- sqrtv * Xp[, -1] + (1 - sqrtv) * Xpmean  
      # Shrink each vector of parameters by the corresponding element in v
      Xp.1 <- Xp[, -1] %*% diag(sqrtv) + X0 %*% ((1 - sqrtv) * Xpmean)
      # Calculate the predicted loss given the shrunk parameter vectors
      Ltilde <- cbind(X0, Xp.1) %*% Bp 
      evpsi[nSamp, l] <- mean(rowMaxs(Ltilde)) #compute EVPSI similar to EVPPI.
    }
  }
  return(evpsi)
}

plotEVPSI <- function(evpsi,
                      parm = "Subset",
                      title = "Expected Value of Partial Sample Information", 
                      txtsize = 12){
  evpsi.df <- data.frame(Parameter = parm, 
                         melt(evpsi,
                              varnames = c("N", "WTP")))
  # Make WTP a factor 
  evpsi.df$WTP <- factor(evpsi.df$WTP)
  # Format WTP labels using $/QALY
  levels(evpsi.df$WTP) <- paste(dollar(as.numeric(colnames(evpsi))), "/QALY", sep = "")
  # Name EVSPI column
  colnames(evpsi.df)[4] <- "EVPSI"
  # Facet EVPSI by WTP
  ggplot(data = evpsi.df, aes(x = N, y = EVPSI)) +
    facet_wrap(~ WTP, ncol = 2) + # , scales = "free_y"
    geom_point() +
    geom_line() +
    ggtitle(title) + 
    scale_x_continuous(breaks = number_ticks(6), labels = comma)+
    scale_y_continuous(breaks = number_ticks(5), labels = comma)+    
    xlab("Sample size (n)") +
    ylab("EVPPI ($)") +
    theme_bw() +
    theme(legend.title = element_text(size = txtsize),
          legend.key = element_rect(colour = "black"),
          legend.text = element_text(size = txtsize),
          title = element_text(face="bold", size=15),
          axis.title.x = element_text(face="bold", size=txtsize),
          axis.title.y = element_text(face="bold", size=txtsize),
          axis.text.y = element_text(size=txtsize),
          axis.text.x = element_text(size=txtsize))
}

COSS <- function(wtp, Strategies, Outcomes, Parameters, tot.pop, 
                 selParams, n0, samp.size, cost.res, parm, ...){
  # Computes COSS for a set of parameters and specific research studies
  # sample sizes
  #
  # Args:
  #   wtp:        Vector of willingness to pay values
  #   Strategies: Vector with strategies' names
  #   Outcomes:   Matrix with outcomes ordered in a way that for each strategy the cost must 
  #               appear first then the effectiveness
  #   Parameters: Matrix with parameters' simulations
  #   tot.pop:    Total population afected by technology
  #   selParams:  Vector of subset of parameters to aggregate into one evpsi measure 
  #               for each WTP at each samp.size
  #   n0:         Vector with effective sample size of parameters
  #   samp.size:  Research study sample sizes
  #   cost.res:   Matirx with cost of research by study sample size (number of rows 
  #               same length as samp.size) and WTP thresholds (number of columns
  #               same as length of wtp)
  #
  # Returns:
  #   coss: Data frame with optimal sample size for different values of WTP
  #
  # Load dependencies 
  require(matrixStats)  # for colMins...etc
  require(reshape2)
  require(dplyr)
  # Generate internal variables
  nStrategy <- length(Strategies) # Number of strategies
  X         <- Parameters # Matrix with input parameters
  nParam    <- ncol(X) # Number of inputs
  nSim      <- nrow(X) # Number of simulations
  nameParam <- colnames(X) # Name of input parameters
  nSample   <- length(samp.size)  # Number of samples
  nWTP      <- length(wtp)  # Number of WTP
  # cost.study <- data.frame(N = samp.size, CS = cost.res)  # Cost of study by sample size
  cost.study <- cost.res
  dimnames(cost.study) <- list(samp.size, wtp)
  
  costInd <- seq(1, 2*nStrategy, by = 2) #vector to index costs
  effInd  <- seq(2, 2*nStrategy, by = 2) #vector to index effectiveness
  evpsi   <- array(0, dim = c(nSample, nWTP),
                   dimnames = list(samp.size, wtp))
  X0 <- matrix(1, nSim, 1) # Intercept for LRM
  Xp <- data.matrix(cbind(X0, X[, selParams])) # Selected vectors
  
  for(l in 1:nWTP){ # l <- 1
    NMB <-  wtp[l]*Outcomes[, effInd]-Outcomes[, costInd] # Effectiveness minus Costs, with vector indexing
    ## Determine the optimal strategy with current info
    dStar <- which.max(colMeans(NMB)) #Compute E[NMB] for all strategies and determine the optimal
    ## Calculate the opportunity loss from choosing dStar
    Loss <- matrix(0, nSim, nStrategy)
    for (s in 1:nStrategy){
      Loss[, s] <- NMB[, s] - NMB[, dStar]
    }
    # Mean of selected vectors
    Xpmean <- colMeans(Xp[, -1])  
    Bp <- solve(t(Xp) %*% Xp) %*% t(Xp) %*% Loss  # Regression coefficient matrix
    for (nSamp in 1:nSample){ # nSamp <- 1
      v      <- samp.size[nSamp] / (samp.size[nSamp] + n0)
      sqrtv  <- sqrt(v)
      Xp.1 <- sqrtv * Xp[, -1] + (1 - sqrtv) * Xpmean  #shrink each vector of parameters by the corresponding element in v
      Xp.1 <- Xp[, -1] %*% diag(sqrtv) + X0 %*% ((1 - sqrtv) * Xpmean)   #shrink each vector of parameters by the corresponding element in v
      Ltilde <- cbind(X0, Xp.1) %*% Bp #calculate the predicted loss given the shrunk parameter vectors
      evpsi[nSamp, l] <- mean(rowMaxs(Ltilde)) #compute EVPSI similar to EVPPI.
    }
  }
  # Use EVSI to compute ENBS
  cost.study.lng <- melt(cost.study, 
                         varnames = c("N", "WTP"), 
                         value.name = "CS")
  enbs <- data.frame(Parameter = parm, # Reshape EVSI into long format and rename as ENBS
                         melt(evpsi,
                              varnames = c("N", "WTP")))  
  colnames(enbs)[4] <- "EVPSI"
  enbs$PopEVPSI <- enbs$EVPSI*tot.pop  # Compute population EVSI
  enbs <- merge(enbs, cost.study.lng, by = c("N", "WTP"))  # Merge with cost of study by sample size 
  enbs <- arrange(enbs, WTP, N)  # Sort by WTP and N
  enbs$ENBS <- enbs$PopEVPSI - enbs$CS  # Compute ENBS
  enbs <- group_by(enbs, WTP)  # Group by WTP
  # OSS: Find at which WTP the ENBS is maximum
  oss <- summarise(enbs,
                   OSS = as.numeric(N[which.max(ENBS)]))
  coss <- data.frame(Parameter = parm, oss)
  return(coss)
}

ENBS <- function(wtp, Strategies, Outcomes, Parameters, tot.pop, 
                 selParams, n0, samp.size, cost.res, parm, ...){
  # Computes COSS for a set of parameters and specific research studies
  # sample sizes
  #
  # Args:
  #   wtp:        Vector of willingness to pay values
  #   Strategies: Vector with strategies' names
  #   Outcomes:   Matrix with outcomes ordered in a way that for each strategy the cost must 
  #               appear first then the effectiveness
  #   Parameters: Matrix with parameters' simulations
  #   tot.pop:    Total population afected by technology
  #   selParams:  Vector of subset of parameters to aggregate into one evpsi measure 
  #               for each WTP at each samp.size
  #   n0:         Vector with effective sample size of parameters
  #   samp.size:  Research study sample sizes
  #   cost.res:   Cost of research by study sample size (same length as samp.size)
  #
  # Returns:
  #   coss: Data frame with optimal sample size for different values of WTP
  #
  # Load dependencies 
  require(matrixStats)  # for colMins...etc
  require(reshape2)
  require(dplyr)
  # Generate internal variables
  nStrategy <- length(Strategies) # Number of strategies
  X         <- Parameters # Matrix with input parameters
  nParam    <- ncol(X) # Number of inputs
  nSim      <- nrow(X) # Number of simulations
  nameParam <- colnames(X) # Name of input parameters
  nSample   <- length(samp.size)  # Number of samples
  nWTP      <- length(wtp)  # Number of WTP
  # cost.study <- data.frame(N = samp.size, CS = cost.res)  # Cost of study by sample size
  cost.study <- cost.res
  dimnames(cost.study) <- list(samp.size, wtp)
  
  costInd <- seq(1, 2*nStrategy, by = 2) #vector to index costs
  effInd  <- seq(2, 2*nStrategy, by = 2) #vector to index effectiveness
  evpsi   <- array(0, dim = c(nSample, nWTP),
                   dimnames = list(samp.size, wtp))
  X0 <- matrix(1, nSim, 1) # Intercept for LRM
  Xp <- data.matrix(cbind(X0, X[, selParams])) # Selected vectors
  
  for(l in 1:nWTP){
    NMB <-  wtp[l]*Outcomes[, effInd]-Outcomes[, costInd] # Effectiveness minus Costs, with vector indexing
    ## Determine the optimal strategy with current info
    dStar <- which.max(colMeans(NMB)) #Compute E[NMB] for all strategies and determine the optimal
    ## Calculate the opportunity loss from choosing dStar
    Loss <- matrix(0, nSim, nStrategy)
    for (s in 1:nStrategy){
      Loss[, s] <- NMB[, s] - NMB[, dStar]
    }
    # Mean of selected vectors
    Xpmean <- colMeans(Xp[, -1])  
    Bp <- solve(t(Xp) %*% Xp) %*% t(Xp) %*% Loss  # Regression coefficient matrix
    for (nSamp in 1:nSample){
      v      <- samp.size[nSamp] / (samp.size[nSamp] + n0)
      sqrtv  <- sqrt(v)
      Xp.1 <- sqrtv * Xp[, -1] + (1 - sqrtv) * Xpmean  #shrink each vector of parameters by the corresponding element in v
      Xp.1 <- Xp[, -1] %*% diag(sqrtv) + X0 %*% ((1 - sqrtv) * Xpmean)   #shrink each vector of parameters by the corresponding element in v
      Ltilde <- cbind(X0, Xp.1) %*% Bp #calculate the predicted loss given the shrunk parameter vectors
      evpsi[nSamp, l] <- mean(rowMaxs(Ltilde)) #compute EVPSI similar to EVPPI.
    }
  }
  # Use EVSI to compute ENBS
  cost.study.lng <- melt(cost.study, 
                         varnames = c("N", "WTP"), 
                         value.name = "CS")
  enbs <- data.frame(Parameter = parm, # Reshape EVSI into long format and rename as ENBS
                     melt(evpsi,
                          varnames = c("N", "WTP")))  
  colnames(enbs)[4] <- "EVPSI"
  enbs$PopEVPSI <- enbs$EVPSI*tot.pop  # Compute population EVSI
  enbs <- merge(enbs, cost.study.lng, by = c("N", "WTP"))  # Merge with cost of study by sample size 
  enbs <- arrange(enbs, WTP, N)  # Sort by WTP and N
  enbs$ENBS <- enbs$PopEVPSI - enbs$CS  # Compute ENBS
  enbs <- group_by(enbs, WTP)  # Group by WTP
  # OSS: Find at which WTP the ENBS is maximum
  oss <- summarise(enbs,
                   OSS = as.numeric(N[which.max(ENBS)]))
  enbs <- summarise(enbs,
                   ENBS = as.numeric(ENBS[which.max(ENBS)]))
  
  
  enbs <- data.frame(Parameter = parm, oss, enbs)
  return(enbs)
}


#### Formatting functions ####
number_ticks <- function(n) {function(limits) pretty(limits, n)} #Function for number of ticks in ggplot
# write a simple function to add footnote
makeFootnote <- function(footnoteText =format(Sys.time(), "%d %b %Y"),size = .7, color = grey(.5)){
  require(grid)
  pushViewport(viewport())
  grid.text(label = footnoteText ,
            x = unit(1,"npc") - unit(2, "mm"),
            y = unit(2, "mm"),
            just = c("right", "bottom"),
            gp = gpar(cex = size, col = color))
  popViewport()
}
# Obtain ggplot default colors to manually add them later on
gg_color_hue <- function(n) { 
  hues = seq(15, 375, length = n+1)
  hcl(h = hues, l = 50, c = 100)[1:n]
}