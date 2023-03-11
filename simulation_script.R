
library(cna)
library(stringr)
library(parallel)

set.seed(677865)

#################
### functions ###
#################


# Function to manipulate prevalence.
introInflation.with.constNoise <- function(datNoisy, groundTruth, inflate.ratio, t, constant = c("sample", "fragmentation"))
{
  datInflate <- vector("list", length(datNoisy))
  if(constant=="sample"){
    for(i in 1:length(t)){
      if(t[[i]] == toupper(t[[i]])) {
        noise.natural <- round(checkNoise(GT=groundTruth[i],x=datNoisy[[i]]),2)
        r <- condition(groundTruth[i], datNoisy[[i]])
        r <- ct2df(r[[1]])
        difference <- apply(r, 1, function(y) abs(y[1]-y[2]))
        compatible <- datNoisy[[i]][which(difference<=noise.natural+0.15),]
        incompatible <- dplyr::setdiff(datNoisy[[i]],compatible)
        select1.compatible <-  subset(compatible, compatible[,which(names(compatible)%in% t[[i]])]>=0.5)
        select1.incompatible <-  subset(incompatible, incompatible[,which(names(incompatible)%in% t[[i]])]>=0.5)
        select1 <- rbind(select1.compatible,select1.incompatible)
        select0.compatible <-  subset(compatible, compatible[,which(names(compatible)%in% t[[i]])]<0.5)
        select0.incompatible <-  subset(incompatible, incompatible[,which(names(incompatible)%in% t[[i]])]<0.5)
        select0 <- rbind(select0.compatible,select0.incompatible)
        datA1.compatible <- some(select1.compatible, round(c*inflate.ratio*1-noise.natural), replace = T)
        datA1.incompatible <- some(select1.incompatible, round(c*inflate.ratio*(noise.natural)), replace = T)
        datA0.compatible <- some(select0.compatible, round(c*(1-inflate.ratio)*1-noise.natural), replace = T)
        datA0.incompatible <- some(select0.incompatible, round(c*(1-inflate.ratio)*(noise.natural)), replace = T)
        datInflate[[i]] <- rbind(datA1.compatible, datA1.incompatible, datA0.compatible, datA0.incompatible) 
        
      } else {
        noise.natural <- round(checkNoise(GT=groundTruth[[i]],x=datNoisy[[i]]),2)
        r <- condition(groundTruth[i], datNoisy[[i]])
        r <- ct2df(r[[1]])
        difference <- apply(r, 1, function(y) abs(y[1]-y[2]))
        compatible <- datNoisy[[i]][which(difference<=noise.natural+0.15),]
        incompatible <- dplyr::setdiff(datNoisy[[i]],compatible)
        select1.compatible <-  subset(compatible, compatible[,which(names(compatible)%in% toupper(t[[i]]))]<0.5)
        select1.incompatible <-  subset(incompatible, incompatible[,which(names(incompatible)%in% toupper(t[[i]]))]<0.5)
        select1 <- rbind(select1.compatible,select1.incompatible)
        select0.compatible <-  subset(compatible, compatible[,which(names(compatible)%in% toupper(t[[i]]))]>=0.5)
        select0.incompatible  <-  subset(incompatible, incompatible[,which(names(incompatible)%in% toupper(t[[i]]))]>=0.5)
        select0 <- rbind(select0.compatible,select0.incompatible)
        datA1.compatible <- some(select1.compatible, round(c*inflate.ratio*1-noise.natural), replace = T)
        datA1.incompatible <- some(select1.incompatible, round(c*inflate.ratio*(noise.natural)), replace = T)
        datA0.compatible <- some(select0.compatible, round(c*(1-inflate.ratio)*1-noise.natural), replace = T)
        datA0.incompatible <- some(select0.incompatible, round(c*(1-inflate.ratio)*(noise.natural)), replace = T)
        datInflate[[i]] <- rbind(datA1.compatible, datA1.incompatible, datA0.compatible, datA0.incompatible) 
        
      }
    }
  }
  if(constant=="fragmentation"){
    for(i in 1:length(t)){
      
      if(t[[i]] == toupper(t[[i]])) {
        noise.natural <- round(checkNoise(GT=groundTruth[i],x=datNoisy[[i]]),2)
        r <- condition(groundTruth[i], datNoisy[[i]])
        r <- ct2df(r[[1]])
        difference <- apply(r, 1, function(y) abs(y[1]-y[2]))
        compatible <- datNoisy[[i]][which(difference<=noise.natural+0.15),]
        incompatible <- dplyr::setdiff(datNoisy[[i]],compatible)
        select1.compatible <-  subset(compatible, compatible[,which(names(compatible)%in% t[[i]])]>=0.5)
        select1.incompatible <-  subset(incompatible, incompatible[,which(names(incompatible)%in% t[[i]])]>=0.5)
        select1 <- rbind(select1.compatible,select1.incompatible)
        select0.compatible <-  subset(compatible, compatible[,which(names(compatible)%in% t[[i]])]<0.5)
        select0.incompatible <-  subset(incompatible, incompatible[,which(names(incompatible)%in% t[[i]])]<0.5)
        select0 <- rbind(select0.compatible,select0.incompatible)
        inflate.natural <- nrow(select1)/nrow(datNoisy[[i]])
        
        if(inflate.ratio>inflate.natural){
          draw <- ((nrow(datNoisy[[i]])*inflate.ratio) - nrow(select1))/(1-inflate.ratio)
          datA1.compatible <- some(select1.compatible, draw*(1-noise.natural), replace = T)
          datA1.incompatible <- some(select1.incompatible, draw*(noise.natural), replace = T)
          datInflate[[i]]  <- rbind(datA1.compatible, datA1.incompatible, datNoisy[[i]]) 
        }else{
          draw <- ((nrow(datNoisy[[i]])*(1-inflate.ratio)) - nrow(select0))/(1-(1-inflate.ratio))
          datA0.compatible <- some(select0.compatible, draw*(1-noise.natural), replace = T)
          datA0.incompatible <- some(select0.incompatible, draw*(noise.natural), replace = T)
          datInflate[[i]]  <- rbind(datA0.compatible,datA0.incompatible, datNoisy[[i]])        
        }
        
      } else {
        noise.natural <- round(checkNoise(GT=groundTruth[[i]],x=datNoisy[[i]]),2)
        r <- condition(groundTruth[i], datNoisy[[i]])
        r <- ct2df(r[[1]])
        difference <- apply(r, 1, function(y) abs(y[1]-y[2]))
        compatible <- datNoisy[[i]][which(difference<=noise.natural+0.15),]
        incompatible <- dplyr::setdiff(datNoisy[[i]],compatible)
        select1.compatible <-  subset(compatible, compatible[,which(names(compatible)%in% toupper(t[[i]]))]<0.5)
        select1.incompatible <-  subset(incompatible, incompatible[,which(names(incompatible)%in% toupper(t[[i]]))]<0.5)
        select1 <- rbind(select1.compatible,select1.incompatible)
        select0.compatible <-  subset(compatible, compatible[,which(names(compatible)%in% toupper(t[[i]]))]>=0.5)
        select0.incompatible  <-  subset(incompatible, incompatible[,which(names(incompatible)%in% toupper(t[[i]]))]>=0.5)
        select0 <- rbind(select0.compatible,select0.incompatible)
        
        inflate.natural <- nrow(select1)/nrow(datNoisy[[i]])
        
        if(inflate.ratio>inflate.natural){
          draw <- ((nrow(datNoisy[[i]])*inflate.ratio) - nrow(select1))/(1-inflate.ratio)
          datA1.compatible <- some(select1.compatible, draw*(1-noise.natural), replace = T)
          datA1.incompatible <- some(select1.incompatible, draw*(noise.natural), replace = T)
          datInflate[[i]]  <- rbind(datA1.compatible, datA1.incompatible, datNoisy[[i]]) 
        }else{
          draw <- ((nrow(datNoisy[[i]])*(1-inflate.ratio)) - nrow(select0))/(1-(1-inflate.ratio))
          datA0.compatible <- some(select0.compatible, draw*(1-noise.natural), replace = T)
          datA0.incompatible <- some(select0.incompatible, draw*(noise.natural), replace = T)
          datInflate[[i]]  <- rbind(datA0.compatible,datA0.incompatible, datNoisy[[i]]) 
          
        }
        
      }
    }
    
  }
  datInflate
}


# Function to check noise (auxiliary for introInflation.with.constNoise).
checkNoise <- function(x, GT){
  r <- condition(GT, x)
  r <- ct2df(r[[1]])
  r <- round(r,2)
  # apply(r, 1, function(y) abs(y[1]-y[2])) %>% table
  out <- mean(apply(r, 1, function(y) abs(y[1]-y[2])))
  out
}


# Function to remove brackets.
remove_par <- function(string){
  gsub("^\\(|\\)$", "", string)
}


# Function to generate dataset, conduct CNA analysis, and process returned models.
run_analysis <- function(seed, setting){
  # Set seed.
  set.seed(seed = seed)
  # Assign settings.
  con <- setting$con_threshold
  cov <- setting$cov_threshold
  factors <- unlist(setting$factors)
  frag_ratio <- setting$frag_ratio
  if (frag_ratio == -99){
    frag_ratio <- sample(seq(0.2,0.5,0.01),1)
  }
  noise_level <- setting$noise_level
  if (noise_level == -99){
    noise_level <- sample(seq(0.05,0.3,0.01),1)
  }
  # Assign prevalence.
  prevalence <- setting$prevalence
  
  # Start generating dataset: get all logically possible variable value combinations.
  fullData <- allCombs(factors) - 1
  # Draw a DGCS (ground truth), with maximum number of conjuncts per conjunction 
  # and disjuncts per disjunction (compl) sampled from 2:4.
  if (outcome != TRUE){
    groundTruth <- randomCsf(fullData, n.asf = 3, compl = 2:4, outcome = outcome)
  }else{
    groundTruth <- randomCsf(fullData, n.asf = 3, compl = 2:4)
  }
  # Generate ideal data for the DGCS (dataset containing all cases compatible 
  # with the DGCS once).
  idealData <- ct2df(selectCases(groundTruth, fullData))
  # Introduce frag_ratio fragmentation.
  if (round(nrow(idealData)*frag_ratio) == 0){
    fragData <- idealData
  }else{
    fragData <- idealData[-sample(1:nrow(idealData), round(nrow(idealData)*frag_ratio)), ]
  }
  # Sample with replacement from the fragmented data.
  if (sample_size != FALSE){
    sampledData <- fragData[sample(1:nrow(fragData), sample_size, replace = TRUE), ]
  }else{
    sampledData <- fragData
  }
  # Get cases incompatible with DGCS (noisy cases)
  incompCases <- dplyr::setdiff(fullData, idealData)
  # Introduce noise_level random noise (cases incompatible with ground truth).
  noisyData <- rbind(incompCases[sample(1:nrow(incompCases),
                                        + round(nrow(sampledData) * noise_level), replace = TRUE), ], sampledData)
  # Reformat DGCS for use as argument in introInflation.with.constNoise().
  theGT <- remove_par(toString(groundTruth))
  # Change prevalence to required level while keeping fragmentation and noise
  # constant.
  x <- introInflation.with.constNoise(datNoisy=list(noisyData), 
                                      groundTruth = theGT, t = list("A"), 
                                      inflate.ratio = prevalence, constant = "fragmentation")[[1]]
  
  # Run CNA on generated dataset.
  nrow_csfs <- 0
  csfs <- tryCatch({
    csfs <- csf(cna(x, outcome = outcome, con = con, cov = cov, maxstep = c(4, 4, 16), rm.dup.factors = F,
                    rm.const.factors = F))
  }, error = function(err)  {
    cat(paste("Error\n Seed:", as.character(seeds[i]),
              "\n Fragmentation ratio:", as.character(frag_ratio),
              "\n Noise level:", as.character(noise_level), "\n"))
    csfs <- list()
  })
  # Return -9 values if no models are returned by CNA (these results will be
  # discarded).
  if (nrow(csfs) < 1){
    result <- rep(-9,14)
    return(result)
  }
  
  # Get sample size of analyzed dataset (as sample size may have been changed by adjusting prevalence
  # using introInflation.with.constNoise).
  N <- nrow(x)
  # Get |Y|.
  Y <- sum(x[csfs$outcome[1]])
  # Get |y|.
  y <- N - Y
  # Get TP, or |X*Y|.
  csfs$TP <- csfs$coverage*Y
  # Get FN, or |x*Y|.
  csfs$FN <- Y - csfs$TP
  # Get |X|.
  csfs$X <- csfs$TP/csfs$consistency
  # Get |x|.
  csfs$x_small <- N - csfs$X
  # Get FP, or |X*y|.
  csfs$FP <- csfs$X - csfs$TP
  # Get TN, or |x*y|.
  csfs$TN <- N - csfs$TP - csfs$FN - csfs$FP
  # Get contrapositive consistency.
  csfs$contr_con <- csfs$TN/y
  # Get contrapositive coverage.
  csfs$contr_cov <- csfs$TN/csfs$x_small
  
  # Get correctness of each model.
  csfs$correct <- is.submodel(csfs$condition, groundTruth)
  # Get DGCS complexity.
  GT_compl <- str_count(lhs(groundTruth),"[A-Za-z]")
  # Get degree of completeness of each model.
  # (Degree of completeness is NaN if there are no correct models.)
  csfs$complx <- csfs$correct*csfs$complexity/GT_compl
  # Get absolute completeness of each model (true if model identical to DGCS, false
  # otherwise).
  # (This was recorded in the experiment but not reported in the paper.)
  csfs$complete <- lapply(csfs$condition, function(x, y) identical.model(x, groundTruth))
  
  # Get low contrapositive models (models not reaching both contrapositive
  # thresholds).
  low_contr_AND <- subset(csfs, contr_con < con | contr_cov < cov)
  # Get correctness of low contrapositive models.
  low_contr_AND_correct <- mean(as.numeric(low_contr_AND$correct))
  # Get degree of completeness of low contrapositive models.
  low_contr_AND_complx <- mean(as.numeric(low_contr_AND$complx))/low_contr_AND_correct
  # Get absolute completeness of low contrapositive models. (Not reported in 
  # paper.)
  low_contr_AND_complete <- mean(as.numeric(low_contr_AND$complete))
  # Get model complexity of low contrapositive models.
  low_contr_AND_nct <- mean(as.numeric(low_contr_AND$complexity))
  
  # Get high contrapositive models. (Models reaching both contrapositive
  # thresholds.)
  high_contr_AND <- subset(csfs, contr_con >= con & contr_cov >= cov)
  # Get correctness of high contrapositive models.
  high_contr_AND_correct <- mean(as.numeric(high_contr_AND$correct))
  # Get degree of completeness of high contrapositive models.
  high_contr_AND_complx <- mean(as.numeric(high_contr_AND$complx))/high_contr_AND_correct
  # Get absolute completeness of high contrapositive models. (Not reported in 
  # paper.)
  high_contr_AND_complete <- mean(as.numeric(high_contr_AND$complete))
  # Get model complexity of high contrapositive models.
  high_contr_AND_nct <- mean(as.numeric(high_contr_AND$complexity))
  
  # Get the ratio of models in the output that are high contrapositive models.
  high_contr_AND_ratio <- nrow(high_contr_AND)/nrow(csfs)
  
  # Get correctness of all returned models (both low and high contrapositive
  # models).
  base_correct <- mean(as.numeric(csfs$correct))
  # Get degree of completeness of all returned models.
  base_complx <- mean(as.numeric(csfs$complx))/base_correct
  # Get absolute completeness of all returned models.
  base_complete <- mean(as.numeric(csfs$complete))
  # Get model complexity of all returned models.
  base_nct <- mean(as.numeric(csfs$complexity))

  # Record the results to return.
  result <- list(low_contr_AND_correct, high_contr_AND_correct, 
                 low_contr_AND_complx, high_contr_AND_complx,
                 low_contr_AND_complete, high_contr_AND_complete,
                 low_contr_AND_nct, high_contr_AND_nct,
                 high_contr_AND_ratio,
                 base_correct, base_complx, base_complete,
                 base_nct, GT_compl
  )
  return(result)
}



############
### main ###
############


# Initialize the dataframe for recording the results of the experiment.
info_df <- data.frame(matrix(ncol =16, nrow=0))

# Set the number of datasets to be generated and analyzed for each combination
# of settings.
n <- 50000

# Sample the seeds for the experiment.
seeds <-  sample(1:1000000, n)

# Set the outcome variable value for the experiment.
outcome <- "A"
# Set the sample size for each dataset to be analyzed by CNA in the experiment.
sample_size <- 100
# Choose the settings for which the experiment is conducted.
settings <- expand.grid(con_threshold = seq(0.7, 0.85, 0.05),
                        cov_threshold = seq(0.7, 0.85, 0.05),
                        noise_level = seq(-99, -99, 0.1),
                        frag_ratio = seq(-99, -99, 0.1),
                        prevalence = seq(0.6, 0.9, 0.05),
                        factors = list(rep(2, 7)))

# Set the number of cores to use for running the experiment.
numCores <- detectCores()

# Run the experiment for each combination of settings.
for (i in 1:nrow(settings)){
  # Apply run_analysis() for each seed for the current combination of settings.
  output <- mclapply(seeds, run_analysis, setting = settings[i,], mc.cores = numCores)
  # Convert the output to a dataframe.
  info_df <- do.call(rbind.data.frame, output)
  # Name the results.
  names(info_df) <- c("low_contr_AND_correct", "high_contr_AND_correct", 
                      "low_contr_AND_complx", "high_contr_AND_complx", 
                      "low_contr_AND_complete", "high_contr_AND_complete", 
                      "low_contr_AND_nct", "high_contr_AND_nct", 
                      "high_contr_AND_ratio", "base_correct", 
                      "base_complx", "base_complete", "base_nct",
                      "GT_complx"
  )
  # Count for how many seeds no CNA models were returned.
  few_csfs <- sum(as.numeric(info_df$base_correct == -9))
  # Discard results for which no CNA models were returned.
  info_df2 <- subset(info_df, base_correct > -9)
  # Count for how many seeds only low contrapositive models were returned.
  all_low_contr_AND <- sum(as.numeric(is.nan(info_df2$high_contr_AND_correct)))
  # Count for how many seeds only high contrapositive models were returned.
  all_high_contr_AND <- sum(as.numeric(is.nan(info_df2$low_contr_AND_correct)))
  # Count for how many seeds complexity results were recorded for low contrapositive models.
  n_low_contr_complx <- sum(as.numeric(!is.nan(info_df2$low_contr_AND_complx)))
  # Count for how many seeds complexity results were recorded for high contrapositive models.
  n_high_contr_complx <- sum(as.numeric(!is.nan(info_df2$high_contr_AND_complx)))
  # Get DGCS complexity of high contrapositive models for each seed.
  GT_compl_h <- info_df2[!is.nan(info_df2$high_contr_AND_correct),]$GT_complx
  # Get DGCS complexity of high contrapositive models averaged over all seeds.
  GT_compl_high_contr <- ifelse(!is.null(GT_compl_h), mean(GT_compl_h), NaN)
  # Get DGCS complexity of low contrapositive models.
  GT_compl_l <- info_df2[!is.nan(info_df2$low_contr_AND_correct),]$GT_complx
  # Get DGCS complexity of low contrapositive models averaged over all seeds.
  GT_compl_low_contr <- ifelse(!is.null(GT_compl_l), mean(GT_compl_l), NaN)
  # Add means over all seeds of recorded results to final dataframe.
  settings[i,names(info_df2)] <- colMeans(info_df2, na.rm = TRUE)
  # Add standard deviations of selected columns to final dataframe.
  settings[i,c("low_corr_sd", "high_corr_sd",
               "low_complx_sd", "high_complx_sd",
               "low_complete_sd", "high_complete_sd",
               "base_corr_sd", "base_complx_sd",
               "base_complete_sd")] <- apply(info_df2[,-c(7,8,9,13,14)], 2, sd, na.rm = TRUE)
  # Add other information to final dataframe.
  settings[i,c("nrow", "few_csfs", "all_high_contr_AND", "all_low_contr_AND",
               "n_low_contr_complx", "n_high_contr_complx")] <- 
    c(nrow(info_df2), few_csfs, all_high_contr_AND, all_low_contr_AND,
      n_low_contr_complx, n_high_contr_complx)
  # Add DGCS complexity of high and low contrapositive models to final dataframe.
  settings[i,c("GT_high_complx", "GT_low_complx")] <-
    c(GT_compl_high_contr, GT_compl_low_contr)
}

# Change factors from vector length of factor (number of variable values) to be
# able to write to csv file.
settings$factors <- lapply(settings$factors, length)
settings <- apply(settings,2,as.numeric)

write.csv(settings, "discr_results/disc_contr_19jan.csv")
