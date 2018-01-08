runner_single_iteration<-function(data, dependence_type, bias_method, bias_params,TargetSampleSize, input_dir, prms, test_type,statistic_type)
{  
  library(pracma)
  source('Compute_Statistic.R')
  source('Compute_Hoeffding_Statistic.R')
  source('Hyperplane_Truncated_Bivariate_Gaussian.R')
  source('estimate_marginals.R')
  source('get_null_distribution.R')
  source('Bootstrap.R')
  source('Get_Dummy_Sample.R')
  source('Tsai_test.R')
  source('Yanivs_estimate_marginals.R')
  source('utilities.R')
  
  #1.################Creating a fixed grid:#######################
  switch(dependence_type,
         'Gaussian'={t<-meshgrid(seq(-2,2,0.1), seq(-2,2,0.1))},
         'Exponential'={t<-meshgrid(seq(0,1.5,0.05), seq(0,1.5,0.05))},
          'Other'={#t<-meshgrid(seq(-2,2,0.1), seq(-2,2,0.1))} )
                    idx_x_min<-which(data[,1]==min(data[,1]))
                    idx_x_max<-which(data[,1]==max(data[,1]))
                    temp_x<-data[-c(idx_x_min,idx_x_max),1]
                    temp_x<-sort(temp_x)
                    
                    idx_y_min<-which(data[,2]==min(data[,2]))
                    idx_y_max<-which(data[,2]==max(data[,2]))
                    temp_y<-data[-c(idx_y_min,idx_y_max),2]
                    temp_y<-sort(temp_y)
                    t<-list(data[-c(idx_x_min,idx_x_max,idx_y_min,idx_y_max),1],data[-c(idx_x_min,idx_x_max,idx_y_min,idx_y_max),2])
                    names(t)<-c('X', 'Y')})
                    #t<-meshgrid(unique(temp_x), unique(temp_y))})#list(data[-c(idx_x_min,idx_x_max,idx_y_min,idx_y_max),1],data[-c(idx_x_min,idx_x_max,idx_y_min,idx_y_max),2])
                    #names(t)<-c('X', 'Y')})#           
  x<-as.vector(t$X)
  y<-as.vector(t$Y)  
  grid_points<-cbind(x,y)
  
  switch (bias_method,
          'Hyperplane_Truncation' = {signs<-bias_params$Hyperplane_prms[1]*x+bias_params$Hyperplane_prms[2]*y+bias_params$Hyperplane_prms[3]
          idx<-which(signs>0)
          grid_points<-grid_points[idx,]})
  #################################################################
  #2. Compute weights matrix W:    
  W=get.biased.sampling.weights(data, dim(data)[1], bias_method, bias_params)
  #################################################################
 switch(test_type,
  'bootstrap'={
  #3. Estimate the marginals (note that at this stage we assume that we know the marginals, that are standard normal):
  bEstimate_Marginals = prms$bEstimate_Marginals
  #marginals<-Yanivs_estimate_marginals(data, 1)
  browser()
  estimate_marginals(data, dependence_type, bias_method, bias_params, bEstimate_Marginals, prms)
  cdfs<-marginals$CDFs
  pdfs<-marginals$PDFs
  
  #marginals<-estimate_marginals(data, dependence_type, bias_method, bias_params, bEstimate_Marginals, prms)
  #pdfs<-marginals$PDF
  #cdfs<-marginals$CDF
  
  #4. Estimate W(x,y)*F_X*FY/normalizing_factor
  null_distribution<-get_null_distribution(data,pdfs,W)
  normalizing_factor<-null_distribution$normalizing_factor
  null_distribution<-null_distribution$null_distribution
  
  #5. Compute the test statistics:
  if(statistic_type == 'HHG' ||statistic_type == 'kendall')
  {
    #1. First compute the statistics based on the original data set:
    True_T=Compute_Statistic(data, statistic_type, grid_points, data, null_distribution)
    
    #True_T=Compute_Hoeffding_Statistic(data, null_distribution)
    #browser()
    #2. Compute statistic for dummy sample:
    statistics_bootstrap=matrix(0,TargetSampleSize,1)
    for(ctr in 1:TargetSampleSize) 
    {
      #Rprof(tmp <- tempfile())
      dummy_sample<-Bootstrap(data, null_distribution)
      #dummy_sample<-Bootstrap_according_to_null(data, null_distribution)
      #Rprof()
      #summaryRprof(tmp)
      
      #browser()
      
      #tryCatch(
      #{
        #browser()
        marginals_bootstrap<-Yanivs_estimate_marginals(dummy_sample, 1)
        cdfs_bootstrap<-marginals_bootstrap$CDFs
        pdfs_bootstrap<-marginals_bootstrap$PDFs
        #browser()
        #marginals_bootstrap<-estimate_marginals(dummy_sample, dependence_type, bias_method, bias_params, bEstimate_Marginals)
        #pdfs_bootstrap<-marginals_bootstrap$PDF
        #cdfs_bootstrap<-marginals_bootstrap$CDF
        
        #3. Compute weights matrix W:    
        W_bootstrap=get.biased.sampling.weights(dummy_sample, dim(dummy_sample)[1], bias_method, bias_params)

        #4. Estimate W(x,y)*F_X*FY/normalizing_factor
        null_distribution_bootstrap<-get_null_distribution(dummy_sample,pdfs_bootstrap,W_bootstrap)
        normalizing_factor_bootstrap<-null_distribution_bootstrap$normalizing_factor
        null_distribution_bootstrap<-null_distribution_bootstrap$null_distribution
        
        #statistics_bootstrap[ctr]<-Compute_Statistic(dummy_sample, statistic_type, grid_points,data, null_distribution)},
        statistics_bootstrap[ctr]<-Compute_Statistic(dummy_sample, statistic_type, grid_points,dummy_sample, null_distribution_bootstrap)
        #browser()
        #Rprof()
        #summaryRprof(tmp)
      #},
       # error= function(e){statistics_bootstrap[ctr]=Inf})
    }
  }
  #browser()
  output<-list(True_T,statistics_bootstrap)
  names(output)<-c('True_T', 'statistics_bootstrap')
  },
  ########################
  'permutations'={
   source('MCMC_Permutations.R')  
   source('prepare_null_mass_table_given_grid.R')
   source('compute_statistic_based_permutaions.R')  
   source('Evaluate_mass_table_by_permutations.R')
   source('Yanivs_estimate_marginals.R')
   source('Evaluate_mass_table_by_cdfs.R')
   
   grid<-create_grid(data, nb_of_points=2)
   number_of_permutations =  TargetSampleSize 
   Permutations=MCMC_Permutations(data,W,number_of_permutations,dim(data)[1])
   temp<-get_expected_value_per_cell(grid, Permutations, data)
   browser()
   grid<-get_observed_value_per_cell(grid, data)
   results<-get_statistic(grid)
   browser()
   
   if(dependence_type=='Gaussian'& file.exists('resources/permutations_test_mass_table.Rdata'))
   {
     load('resources/permutations_test_mass_table.Rdata')
   }else
   {
    if(dependence_type=='Gaussian')
    {browser()
     mass_table<-prepare_null_mass_table_given_grid(bias_params, grid_points)
     dir.create('resources',showWarnings = FALSE)
     save(mass_table, file='resources/permutations_test_mass_table.Rdata')
    }
    
    #Here we do not assume we know the marginals in advance
    if(dependence_type=='Other')
    {
      marginals<-Yanivs_estimate_marginals(data, 1)
      cdfs<-marginals$CDFs
      pdfs<-marginals$PDFs
      null_distribution<-get_null_distribution(data,pdfs,W)
      #browser()
      #Evaluate the quarters expectations table, for each point in the grid:
      grid_permutation<-MCMC_Permutations(data,W,1,dim(data)[1])
      grid_points<-cbind(data[,1], data[grid_permutation,2])
      
      idx_x_min<-which(grid_points[,1]==min(grid_points[,1]))
      idx_x_max<-which(grid_points[,1]==max(grid_points[,1]))      
      idx_y_min<-which(grid_points[,2]==min(grid_points[,2]))
      idx_y_max<-which(grid_points[,2]==max(grid_points[,2]))
      grid_points<-grid_points[-c(idx_x_min,idx_x_max, idx_y_min, idx_y_max),]
      #######################################################
      #expectations_table<-Evaluate_mass_table_by_permutations(data, Permutations, number_of_permutations, grid_points)
      #browser()
      #expectations_table<-Evaluate_mass_table_by_cdfs(data, grid_points, marginal_cdfs)
      #browser()
      #Compute the "true" statistics value:
      #True_T=compute_statistic_based_permutaions(data, statistic_type, dim(data)[1], bias_method, dependence_type, grid_points, expectations_table)
      True_T<-Compute_Statistic(data, statistic_type, grid_points, data, null_distribution$null_distribution)
      #Permutations_2=MCMC_Permutations(data,W,number_of_permutations,dim(data)[1])
      #Compute the statistics value for each permutation:
      statistics_permutations=matrix(0,TargetSampleSize,1)
      statistics_permutations_cor=matrix(0,TargetSampleSize,1)
      #browser()
      for(ctr in 1:TargetSampleSize) 
      {
        #statistics_permutations[ctr]=compute_statistic_based_permutaions(cbind(data[,1],data[Permutations[,ctr],2]), statistic_type, dim(data)[1],bias_method,dependence_type, grid_points, expectations_table)
        #statistics_permutations[ctr]=(1/(number_of_permutations-1))*compute_statistic_based_permutaions(cbind(data[,1],data[Permutations[,ctr],2]), statistic_type, dim(data)[1],bias_method,dependence_type, grid_points, expectations_table)
        statistics_permutations[ctr]<-Compute_Statistic(cbind(data[,1],data[Permutations[,ctr],2]), statistic_type, grid_points, data, null_distribution$null_distribution)
          
      }
    }}
   output<-list(True_T,statistics_permutations)
   names(output)<-c('True_T', 'statistics_permutations')})
    
  return(output)
}









