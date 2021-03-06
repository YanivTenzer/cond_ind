# Simulate data with different type of dependency 
# n - # of observations
# dependence_type - what dependency to simulate (default is linear)
# dependence_params - parameters for dependency (default is making X,Y independent)
# Important! need to add here simulations with biased sampling, which depend on a function W 

simulate_data <- function(n, dependence_type, prms)
{
  switch(dependence_type,   
         'Gaussian' ={library('mvtnorm')
                    mu=c(0,0)
                    rho = prms$rho
                    Sigma=matrix(c(1, rho, rho,1),2,2)
                    data<-rmvnorm(n, mu, Sigma)
                           },
         'Gumbel'={library('copula')
                   gumbel.cop <- gumbelCopula(0)
                   u<- rCopula(n, gumbel.cop)
                   data<-cbind(qnorm(u[,1]),qnorm(u[,2]))},
         'Clayton'={library('copula')
                    clayton.cop <- claytonCopula(0.0001)
                   u<- rCopula(n, clayton.cop)
                   data<-cbind(qnorm(u[,1]),qnorm(u[,2]))},
         'mixture'={
                    library('copula')
                    clayton.cop <- claytonCopula(1)
                    u_1<- rCopula(n, clayton.cop)
                    clayton.cop <- claytonCopula(1)
                    u_2<- rCopula(n, clayton.cop)
                    data<-rbind(u_1,u_2)
                    #library('mvtnorm')
                    #mu=c(0,0)
                    #rho = 0.95
                    #Sigma_1=matrix(c(1, rho, rho,1),2,2)
                    #Sigma_2=matrix(c(1, -rho, -rho,1),2,2)
                    #data<-rbind(rmvnorm(floor(n/2), mu, Sigma_1),rmvnorm(floor(n/2), mu, Sigma_2))
                    },
                    
        'Linear' = cbind(x, y=x+ prms$noise*rnorm(n)),
        'Parabolic' = {x=rnorm(n) 
                       data<-cbind(x, y=4*(x-.5)^2)+rnorm(n)},
        'Sin' =       {x=rnorm(n) 
                       data<-cbind(x, y=sin(abs(x))+1*rnorm(n))},
        'Cubic' = cbind(x, 128*(x-1/3)^3-48*(x-1/3)^3-12*(x-1/3)+10*  prms$noise*rnorm(n)),
        'Exponential'={data<-cbind(0.6-rexp(n,prms$lambda),rexp(n,prms$lambda))}) 
  
  
  return (data)  
}
##########################
Create_Bias <- function(Data,biased_method, bias_params)
{
  n=dim(Data)[1]
  if(biased_method == 'SOFT_CENS') 
  {
    th = bias_params$L2_th
    L2_Norms<-apply(Data,1,function(x) sqrt(x[1]^2+x[2]^2) )
    Probs<-min(1,1/L2_Norms)
    Toss_Coin = runif(n,0,1)
    truncated_data<-Data[Toss_Coin<Probs,]
  }  
  
  if(biased_method== 'Hyperplane_Truncation')
  {
    Hyperplane_prms = bias_params$Hyperplane_prms
    M<-matrix(0,n,1)
    M=Data[,1]*Hyperplane_prms[1]+ Data[,2]*Hyperplane_prms[2]+ Hyperplane_prms[3]
    Signs=sign(M)
    truncated_data<-Data[Signs>0,]
  }
  
  if(biased_method== 'non_truncation')
  {
    w_2<-apply(Data,1, function(x) exp(-sum(abs(x))))
    u<-runif(dim(Data)[1])
    bReject<-ifelse(u*1>w_2, 1, 0)
    truncated_data<-Data[!bReject,]
  }
  
  return (truncated_data);
}
###############################
# Compute the N*N matrix of samplig weights (add more methods in the future )
get.biased.sampling.weights <- function(Data, N, biased_method, bias_params)
{
  W = matrix(0,N,N);
  if(biased_method == 'SOFT_CENS')
  {
    th=bias_params$L2_th
    for (i in 1:N )
    {
      Temp=cbind(rep(Data[i,1],N), Data[,2])
      L2_norm = apply(Temp,1,function(x) sqrt(x[1]^2+x[2]^2))
      W[i,which(L2_norm<th)]=1
      Id<-which(L2_norm>th)
      W[i,Id] = 1/L2_norm[Id]
    }
  }
  #####
  if(biased_method == 'Hyperplane_Truncation')
  {
    Hyperplane_prms = bias_params$Hyperplane_prms
    for (i in 1:N )
    {
      Temp=cbind(rep(Data[i,1],N), Data[,2])
      S = apply(Temp,1,function(x) (Hyperplane_prms[1]*x[1]+Hyperplane_prms[2]*x[2]+Hyperplane_prms[3]))
      W[i,which(sign(S)>0)]=1
    }
  }

  return (W); 
}
