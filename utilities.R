create_grid<-function(data, nb_of_points)
{
  idx_x_min<-which(data[,1]==min(data[,1]))
  idx_x_max<-which(data[,1]==max(data[,1]))
  
  idx_y_min<-which(data[,2]==min(data[,2]))
  idx_y_max<-which(data[,2]==max(data[,2]))
  
  all_points<-cbind(data[-c(idx_x_min,idx_x_max,idx_y_min,idx_y_max),1],data[-c(idx_x_min,idx_x_max,idx_y_min,idx_y_max),2])
  nb_all_points<-dim(all_points)[1]
  idx_points<-sample(nb_all_points, nb_of_points, replace=FALSE)
  points<-data[idx_points,]
  
  x_segments<-c(-Inf, sort(points[,1], decreasing=FALSE), Inf)
  y_segments<-c(-Inf, sort(points[,2], decreasing=FALSE), Inf)
  grid<-meshgrid(x_segments, y_segments)
  
  cells<-list()
  count = 0 
  for(i in seq(2, dim(grid$X)[1],1))
    for(j in seq(2, dim(grid$Y)[1],1))
    {
      cell<-{}
      cell$down_left_point<-c(grid$X[1, i-1], grid$Y[j-1, 1])
      cell$up_right_point<-c(grid$X[1, i], grid$Y[j, 1])
      #cell$expected<-0
      count=count+1
      cells[[count]]<-cell
    }
  return(cells)
}
###############################################################
get_expected_value_per_cell<-function(cells, permutations, data)
{
  
  expected_value_per_cell<-matrix(0, length(cells),1)
  for(i in seq(1, dim(permutations)[2],1))
  {
    browser()
    permuted_data<-cbind(data[,1], data[permutations[,i],2])
    for(j in seq(1, length(cells),1))
    {
      x<-cells[[j]]
      expected_value_per_cell[j]=expected_value_per_cell[j] + length(which(permuted_data[,1]>x$down_left_point[1]&permuted_data[,1]<=x$up_right_point[1] & permuted_data[,2]>x$down_left_point[2]&permuted_data[,2]<=x$up_right_point[2]))
    }
    
  }
  expected_value_per_cell<- expected_value_per_cell/(dim(permutations)[2])
  num_of_non_zero_cells = length(which(expected_value_per_cell!=0))
  out_cells<-list()
  for(j in seq(1, length(cells),1))
  {
    temp_cell<-{}
    temp_cell$down_left_point<-cells[[j]]$down_left_point
    temp_cell$up_right_point<-cells[[j]]$up_right_point
    temp_cell$expected_mass<-expected_value_per_cell[j]
    out_cells[[j]]<-temp_cell
  }
  outputs<-list(out_cells,num_of_non_zero_cells)
  names(outputs)<-c('cells', 'num_of_non_zero_cells')
  return(outputs)
}
####################################################
get_observed_value_per_cell<-function(cells, data)
{
  out_cells<-list()
  for(j in seq(1, length(cells),1))
  {
    x<-cells[[j]]
    temp_cell<-{}
    temp_cell$down_left_point<-cells[[j]]$down_left_point
    temp_cell$up_right_point<-cells[[j]]$up_right_point
    temp_cell$expected_mass<-cells[[j]]$expected_mass
    temp_cell$observed<-length(which(data[,1]>x$down_left_point[1]&data[,1]<=x$up_right_point[1] & data[,2]>x$down_left_point[2]&data[,2]<=x$up_right_point[2]))
    out_cells[[j]]<-temp_cell
  }
  return(out_cells)
}
##################################################
get_statistic<-function(cells)
{
  T = 0
  epsilon = 10^(-6)
  for(j in seq(1, length(cells),1))
  {
    x<-cells[[j]]
    T=T+(x$observed-x$expected_mass)^2/(x$expected_mass+epsilon)
  }
  outputs<-c(T,num_of_non_zero_cells)
  return(outputs)
}