#The RM_SMOTE function
#an example how to call the function
#notice that just one of the "dup_size" or "eIR" should be passed to the function
#the class variable should be sent as the factor variable
#example 1
#   oversampled_dataset <- RM_SMOTE(dt = dataset, target = '1', dup_size = 3, k = 5, thre = 0.025, alph = 0.5, weight_f = 1)
#example 2
#   oversampled_dataset <- RM_SMOTE(dt = dataset, target = '1', eIR = 1, k = 5, thre = 0.025, alph = 0.5, weight_f = 1)
#-------------------------------------------------------------------------------------------------------------------------------#

#the needed libraries
library(stats)
library(rrcov)
library(heplots)
library(meanShiftR)
library(DescTools)

# ===================================================
#calculating the weights for each minority observation
# ---------------------------------------------------
weighting <- function(xxx, thre, alph = 0.5, weight_func = 1)
  # INPUTS:
  # xxx are the Minority class observations;
  # thre is the parameter related to the chi-square distribution to ignore the outlines
  # alph is the MCD parameter which suppose the alpha amount of observations to robustly estimate the mean and cov matrix
  # and weight_func specifies which weighting function should be used (1 for outliers 0 weights, 2 for 1/MD_square, 3 MD_threshold/MD_sqiare)
  # OUTPUTS:
  # The result of the function is a dataframe consisting the original one and 3 new columns consisting
  # MD (Mahalanobis distance to estimated robust mean), weights, and probability of choosing each obs.
{
  #calculating the threshold from the chi-sqaure distribution with specified thre and degree of freedom
  #degree of freedom is just the number of explanatory variables (features) in your dataset
  degree_of_freedom <- ncol(xxx)-1
  MD_threshold <- qchisq((1-thre), df = degree_of_freedom)
  
  #calculating the Mean and Covariance using MCD approach
  temp2 <- rrcov::CovMcd(xxx[,1:degree_of_freedom], alpha=alph)
  mean_minority_MCD <- rrcov::getCenter(temp2)
  sigma_minority_MCD <- rrcov::getCov(temp2)
  
  #calculating the Mahalonobis Distance for all Minority samples using the robustly estimated Mean and Covariance (MCD) 
  xxx$MD <- heplots::Mahalanobis(xxx[,1:degree_of_freedom], mean_minority_MCD, sigma_minority_MCD, method = "mcd")
  
  #calculate the related weights
  #first assign weight 1 to all
  xxx$weights <- 1
  #modifying the weights of the samples farther than the threshold
  for(i in c(which(xxx$MD > MD_threshold)))
  {
    if(weight_func == 1)  {
      xxx[i,'weights'] = 0
    } else if(weight_func == 2) {
      xxx[i,'weights'] = 1/xxx[i,'MD']
    }else if(weight_func == 3)  {
      xxx[i,'weights'] = MD_threshold/xxx[i,'MD']
    }
  }
  
  #calculating the probability of selecting specific observation from minority class
  xxx$prob <- (xxx$weights)/(sum(xxx$weights)) 
  
  #returning the new dataframe with the calculated MD, weights, and probabilities
  return(xxx)
}


# =================================================== this version ignores all the observations with weight 0
# Obtain a set of synthetic observations by MAHA-SMOTE for a set of minority observations.
# ---------------------------------------------------
RM_SMOTE <- function(dt, target, dup_size = 0, eIR = 0, k, thre = 0.01, alph = 0.75, weight_f = 1)
  # INPUTS:
  # dt are the whole observations, the class variable should be existed as 'class';
  # target is the value of class variable that shows the minority class;
  # notice that just one of the "dup_size" or "r" should be passed to the function
  # dup_size is the number of times that over-sampling should increase the number of Minority class observations; 
  #     so if you have 100 Minority obs. and dup_size is 3 at the final oversampled dataset you will have 400 Minority obs. 
  # r is the indicator to say we want the completely balanced set ratio that we want to reach from over-sampling;
  # k is the number of nearest neighours (by MAHALANOBIS distance metric)
  # thre which indicates that the threshold to find the outliers using the chi-square distribution
  #     to use for the generation of synthetic observations;
  # alph is the MCD parameter which suppose the alpha amount of observations to robustly estimate the mean and cov matrix
  # weight_f specifies which weighting function should be used 
  #     (1 for outliers 0 weights, 2 for 1/MD_square, 3 MD_threshold/MD_square)
  # OUTPUTS:
  # The result of the function is a (((N/100)*minority_size)-minority_size) set of generated
  # observation with Minority class values on the target
{
  #reset the Indices of the dataframe
  rownames(dt) <- NULL
  
  #subseting to reach the minority samples
  xxx <- subset(dt, class == target)
  
  ####call the Weighting function####
  #calculating the weights and probability of choosing each minority obs. using specified weighting function
  xxx <- weighting(xxx, thre = thre, alph = alph, weight_func = weight_f)
  rownames(xxx) <- NULL
  
  ####KNN with Mahalanobis Distance####
  #finding the number of variables
  num_col <- ncol(dt)-1
  #finding the Indices of all k-nearest neighbors for all minority samples based on the Mahalanobis Distance
  knn_all_neighbors <- knn_meanShift(xxx[,1:num_col], xxx[,1:num_col], k = min(k+1, NROW(xxx)))
  
  ####Resampling step####
  num_obs <- nrow(dt)
  num_min <- nrow(subset(dt, class == target))
  num_maj <- nrow(subset(dt, class != target))
  #finding the proper N if it is not indicated
  if (dup_size != 0)    {N1 <- round(dup_size*num_min)
  }else if (eIR == 1)      {N1 <- (num_maj-num_min)
  }else if (eIR > 1)       {N1 <- round(num_maj/eIR) - num_min}
  
  #extra steps to ignore the outliers in the resampling step (this would be considered just in case of weight_f equals to 1)
  xxx_1 <- xxx
  xxx <- subset(xxx, weights != 0)
  
  n_min_obs <- nrow(xxx)
  #generating a random sample with regard to estimated probabilities as the first parents of RM_SMOTE
  idx <- sample(1:n_min_obs, N1, replace = TRUE, prob = xxx$prob)
  
  #define an empty data frame to keep new samples 
  new_samples <- data.frame()
  
  #generate synthetic observations
  for (j in c(1:length(idx)))
  {
    #select the first parent as the instances of the random sequence generated in previous step
    parent1 <- xxx[idx[j],1:num_col]
    #generating a random number between 1 to k to choose randomly the second parent from the knn-neighbors of first parent
    temp1 <- sample(1:k, 1)
    #extracting the second parent
    parent2 <- xxx[knn_all_neighbors$neighbors[idx[j],temp1], 1:num_col]
    #generating a random number between 0 and 1 as the coefficents of parents to build the new sample 
    coef2 <- runif(1)
    coef1 <- 1 - coef2
    child <- (coef1 * parent1)+(coef2 * parent2)
    new_samples <- rbind(new_samples, child)
  }
  #assign the label to the synthetic observations
  new_samples$class <- target
  #combining the synthetic observations and original set and
  output_set <- rbind(dt, new_samples)
  row.names(output_set) <- NULL
  
  ####return the synthetic observations attached to orignal dataset as output of the function####
  return(output_set)
}


#-------------------------------------------------------------------------------------------------------------------------------#
####an example with the a KEEL benchmark dataset####
#-------------------------------------------------------------------------------------------------------------------------------#
#read the dataset from your local machine
dt <- read.csv("datasets/ecoli3.csv", header = FALSE)
#assign the name of 'class' to class variable of the dataset
names(dt)[ncol(dt)] <- "class"
dt <- sapply( dt, as.numeric )
dt <- as.data.frame(dt)
#make the class variable as factor variable
dt$class <- as.factor(dt$class)
#drop the binary variables from ecoli3 dataset
dt <- subset(dt, select = -c(3, 4))

#perfrom oversampling
balanced_set <- RM_SMOTE(dt = dt, target = '1', dup_size = 3, k = 5, thre = 0.025, alph = 0.5, weight_f = 1)
balanced_set1 <- RM_SMOTE(dt = dt, target = '1', eIR = 1, k = 5, thre = 0.025, alph = 0.5, weight_f = 1)
#-------------------------------------------------------------------------------------------------------------------------------#
