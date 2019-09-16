############################################
##### This code is used to estimate the coefficients with
##### 1. pooled data without one site "jmdc"
##### 2. pooled data without two sites "jmdc" and "mdcr" at the same time


# function to get pooled estimator
run_pooled <- function(Yall, Xall){
  nam = c("Intercept", colnames(Xall))
  
  fitall = summary(glm(Yall~Xall, family = "binomial"(link = "logit")))
  betaall = fitall$coefficients[,1]
  v_betaall = fitall$cov.scaled
  
  names(betaall) = nam
  colnames(v_betaall) =nam
  rownames(v_betaall) =nam
  
  return(list(beta = betaall,
              var_matrix = v_betaall))
}

######################################
# PART 1.. Pooled estimator (use the pooled data without "jmdc" site) #
# pooled estimates using CCAE + MDCD + MDCR + Optum #
#####################################
Yall_1 = # exclude "jmdc"
Xall_1 = # exclude "jmdc"

# results
result_1 = run_pooled(Yall_1, Xall_1) 

######################################
# PART 2.. Pooled estimator (use the pooled data without "jmdc" and "mdcr" site) #
# pooled estimates using CCAE + MDCD + Optum #
#####################################
Yall_2 = # exclude "jmdc" and "mdcr"
Xall_2 = # exclude "jmdc" and "mdcr"

# results
result_2 = run_pooled(Yall_2, Xall_2)
  
  
