cat("\014")
path = 'D:/cond_ind_local'
setwd(path)
source('Simulate _Data.R')
source('runner_single_iteration.R')
library(foreach)
library(doSNOW)
#library(parallel)
library(doParallel)
library(gdata)

cores=detectCores()
cl<-makeCluster(cores[1]-1) #not to overload your computer registerDoSNOW(cl)
registerDoParallel(cl) 

test_type = 'permutations'#permutations' #'bootstrap' 
statistic_type = 'HHG'#kendall
str<-paste(test_type, '_', statistic_type, sep='')
switch(str,
       'bootstrap_HHG'={output_dir<-'results_bootstrap'},
       'bootstrap_kendall'={output_dir<-'results_bootstrap_kendall'},
       'permutations_HHG'={output_dir<-'results_permutations'})
if(!dir.exists(output_dir))
{
  dir.create(output_dir,showWarnings = FALSE)
}

#The size of the original sample (i.e., before truncation)
sample_size=3000
dependence_type='Gaussian'#Other'#'Gaussian' #Exponential
Hyperplane_prms<-c(-1,1,0)
bias_params<-list(Hyperplane_prms)
names(bias_params)<-c("Hyperplane_prms")
bias_method = 'non_truncation'#'Hyperplane_Truncation'  
bEstimate_Marginals = 1
num_of_iterations = 100
#How many permutations to use/how many bootstrap samples to use:
num_of_repetitions_per_iterations = 100

for(iter in c(6))#seq(1,9,1))
{
  switch(dependence_type,
         'Gaussian'={rho =0.1*iter
                    prms<-list(rho, bEstimate_Marginals)
                    names(prms)<-c("rho", 'bEstimate_Marginals')},
         'Exponential'={lambda=20
                        prms<-list(lambda, bEstimate_Marginals)
                        names(marginal_prms)<-c('lambda', 'bEstimate_Marginals')})       
  
  #results_table<-matrix(0, num_of_iterations, 1+num_of_repetitions_per_iterations)
  #The main loop:
  results_table<- foreach(i=1:num_of_iterations, .combine=rbind) %dopar%{ 
                  data<-simulate_data(sample_size, dependence_type, prms) 
                  truncated_data<-Create_Bias(data, bias_method, bias_params)
                  results<-runner_single_iteration(truncated_data, 'Other', bias_method, bias_params,num_of_repetitions_per_iterations, 
                                                    path, prms, test_type,statistic_type)
   switch(test_type,
   'bootstrap'={temp<-c(results$True_T, results$statistics_bootstrap)},
   'permutations'={temp<-c(results$True_T, results$statistics_permutations)})
    
    temp
  }
  
  str<-paste(test_type, '_', statistic_type, sep='')
  switch(str,
        'bootstrap_HHG'={outfile<-paste(output_dir, '/Bootstrap_results_rho_', toString(iter), '.Rdata', sep='')},
        'bootstrap_kendall'={outfile<-paste(output_dir, '/Bootstrap_kendall_results_rho_', toString(iter), '.Rdata', sep='')},  
        'permutations_HHG'={outfile<-paste(output_dir, '/Permutations_results_rho_', toString(iter), '.Rdata', sep='')})
  
  Pvalues<-apply(results_table, 1, function(x) length(which(x[-1]>x[1]))/length(x))
  save(Pvalues,results_table, file=outfile)
}
###############################################
stopCluster(cl)
###############################################
source('runner_single_iteration.R')
num_of_iterations = 200
sample_size=10000
statistic_type = 'HHG'
test_type = 'permutations'

cores=detectCores()
cl<-makeCluster(cores[1]-1) #not to overload your computer registerDoSNOW(cl)
registerDoParallel(cl) 

switch(statistic_type,
        'HHG'={str<-paste(test_type, '_', statistic_type, '_unknown_margins',sep='')},
       'kendall'={str<-paste(test_type, '_', statistic_type, sep='')})

switch(str,
       'permutations_HHG_unknown_margins'={output_dir<-'results_permutations_HHG_unknown_margins'},
       'permutations_kendall'={output_dir<-'results_bootstrap_kendall'})
if(!dir.exists(output_dir))
{
  dir.create(output_dir,showWarnings = FALSE)
}

for(iter in seq(1,9,1))
{
  dependence_type='Gaussian'
  switch(dependence_type,
       'Gaussian'={rho =0.1*iter
                   prms<-list(rho, bEstimate_Marginals)
                   names(prms)<-c("rho", 'bEstimate_Marginals')},
       'Exponential'={lambda=20
                      prms<-list(lambda, bEstimate_Marginals)
                      names(marginal_prms)<-c('lambda', 'bEstimate_Marginals')})       

  results_table<- foreach(i=1:num_of_iterations, .combine=rbind) %dopar%{ 
                  data<-simulate_data(sample_size, dependence_type, prms) 
                  truncated_data<-Create_Bias(data, bias_method, bias_params)
                  results<-runner_single_iteration(truncated_data, 'Other', bias_method, bias_params,num_of_repetitions_per_iterations, 
                                                    path, prms, test_type,statistic_type)
                  switch(test_type,
                         'bootstrap'={temp<-c(results$True_T, results$statistics_bootstrap)},
                         'permutations'={temp<-c(results$True_T, results$statistics_permutations)})
                  
                  temp }
  
  str<-paste(test_type, '_', statistic_type, sep='')
  switch(str,
         'permutations_HHG'={outfile<-paste(output_dir, '/Bootstrap_results_rho_', toString(iter), '.Rdata', sep='')},
         'permutations_kendall'={outfile<-paste(output_dir, '/Bootstrap_kendall_results_rho_', toString(iter), '.Rdata', sep='')})
  
  Pvalues<-apply(results_table, 1, function(x) length(which(x[-1]>x[1]))/length(x))
  save(Pvalues,results_table, file=outfile)
}