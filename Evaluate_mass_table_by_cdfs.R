#This function get a list of points(grid_points) and compute Expect[Q_i(p_j)] for 1<=i<=4, and all j
#For more detalis see Or's comment in this paper
Evaluate_mass_table_by_cdfs<-function(data, grid_points, marginal_cdfs)
{
  n<-dim(data)[1]
  t<-meshgrid(unique(data[,1]), unique(data[,2]))
  x<-as.vector(t$X)
  y<-as.vector(t$Y)  
  temp_grid<-cbind(x,y)
  mass_table<-matrix(0,dim(grid_points)[1],4)
  
  temp_x<-seq(-5, 5, 0.01)
  temp_y<-seq(-5, 5, 0.01)
  t<-meshgrid(temp_x, temp_y)
  l<-length(temp_x)
  f_x<-(marginal_cdfs$F_1(temp_x[-1])+marginal_cdfs$F_2(temp_x[-1]))/2-(marginal_cdfs$F_1(temp_x[-l])+marginal_cdfs$F_2(temp_x[-l]))
  f_x<-c(f_x, f_x[l])
  f_y<-(marginal_cdfs$F_2(temp_y[-1])+marginal_cdfs$F_1(temp_y[-1]))/2-(marginal_cdfs$F_2(temp_y[-l])+marginal_cdfs$F_1(temp_y[-l]))
  f_y<-c(f_y, f_y[l])
  t_2<-meshgrid(f_x, f_y)
  
  for(i in seq(1, dim(grid_points)[1],1))
  {
    p<-grid_points[i,]
    idx_q1<-which(t[,1]>p[1] & t[,2]>p[2])
    temp<-t[idx_q1,]
    mass<-t_2[,]
    #x<-temp[1]
    #F_x<-(marginal_cdfs$F_1(x)+marginal_cdfs$F_2(x))/2
    #y<-temp[2]
    #F_y<-(marginal_cdfs$F_1(y)+marginal_cdfs$F_2(y))/2
    
    P1<-n*(1-F_x)*(1-F_y)
    P2<-n*(1-F_x)*F_y
    P3<-n*F_x*F_y
    P4<-n*F_x*(1-F_y)
    mass_table[i,]<-c(P1,P2,P3,P4)
  } 
  
  mass_table<-mass_table+0.000001
  return(mass_table)
}
