#-----------------------------------------------------------------------------------
# script to include the batchadjustment into the cross validation
#                                                                     
# Martin Sill
# m.sill@dkfz.de                                                                  
# 
# 2018-03-14 UTC
#------------------------------------------------------------------------------------                   
#options(max.print = 1000)
#options(stringsAsFactors = FALSE)
#options(scipen = 999)


batchadjust <- function(betas,batchall,fold){
  # recalculate betas, illumina like
  betas_out <- betas[,fold$train]
  betas_out <- as.data.frame(t(betas_out))

  # illumina-like beta values
  betas.test_out <- betas[,fold$test]
  betas.test_out <- as.data.frame(t(betas.test_out))
  return(list(betas.train=betas_out,betas.test=betas.test_out))
}
