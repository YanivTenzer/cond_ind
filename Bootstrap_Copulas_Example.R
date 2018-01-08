cat("\014")
path = 'D:/cond_ind_local'
setwd(path)
source('Simulate _Data.R')
source('runner_single_iteration.R')
library(foreach)
library(doSNOW)
library(parallel)
library(doParallel)
library(gdata)
require(gdata)

cores=detectCores()
cl<-makeCluster(cores[1]-1) #not to overload your computer registerDoSNOW(cl)
registerDoParallel(cl) 

test_type = 'bootstrap'
statistic_type = 'HHG'
str<-paste(test_type, '_', statistic_type, sep='')
switch(str,
       'bootstrap_HHG'={output_dir<-'results_bootstrap'},
       'permutations_HHG'={output_dir<-'results_permutations'})
if(!dir.exists(output_dir))
{
  dir.create(output_dir,showWarnings = FALSE)
}

#The size of the original sample (i.e., before truncation)
sample_size=500
dependence_type='mixture'
Hyperplane_prms<-c(-1,1,0)
bias_params<-list(Hyperplane_prms)
names(bias_params)<-c("Hyperplane_prms")
bias_method = 'Hyperplane_Truncation'  
bEstimate_Marginals = 1
num_of_iterations = 100
#How many permutations to use/how many bootstrap samples to use:
num_of_repetitions_per_iterations = 100


results_table<- foreach(i=1:num_of_iterations, .combine=rbind) %dopar%{ 
data<-simulate_data(sample_size, dependence_type, NULL) 
truncated_data<-Create_Bias(data, bias_method, bias_params)
results<-runner_single_iteration(truncated_data, 'Other', bias_method, bias_params,num_of_repetitions_per_iterations, path, NULL, test_type,statistic_type)

switch(test_type,
           'bootstrap'={temp<-c(results$True_T, results$statistics_bootstrap)},
           'permutations'={temp<-c(results$True_T, results$statistics_permutations)})
    
    temp
}

str<-paste(test_type, '_', statistic_type, dependence_type,sep='')
switch(str,
       'bootstrap_HHG'={outfile<-paste(output_dir, '/Bootstrap_results_', dependence_type, '.Rdata', sep='')},  
       'permutations_HHG'={outfile<-paste(output_dir, '/Permutations_results_', dependence_type, '.Rdata', sep='')})

Pvalues<-apply(results_table, 1, function(x) length(which(x[-1]>x[1]))/length(x))
save(Pvalues,results_table, file=outfile)
  