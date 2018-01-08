cat("\014")
rm(list=ls())
path = 'D:/cond_ind_local'
setwd(path)
source('MCMC_Permutations.R')  
########################
#This is a toy example, simulating bivariate vector (x,y) \in {0,1,2}x{0,1,2}):
#1. Define P(X,Y):
probability<-matrix(0,3,3)
probability[1,1]=1/6
probability[2,2]=1/6
probability[3,3]=1/6
probability[1,2]=1/12
probability[1,3]=1/12
probability[2,1]=1/12
probability[2,3]=1/12
probability[3,1]=1/12
probability[3,2]=1/12
#2. W*P(x,y):
biased_probability<-matrix(0,3,3)
#w(0,0)=w(0,1)=w(0,2)=w(1,2)=w(2,1)=0:
S<-sum(colSums(probability))-(sum(probability[,3])+sum(probability[3,]))+probability[3,3]
#Normalize the resulting distribution:
biased_probability[1:2,1:2]<-probability[1:2,1:2]/S
#biased_probability[2,3]=0
#biased_probability[3,2]=0
#3. Sample from  W*P(x,y):
n=20000
rand_idx<-runif(n)
idx<-matrix(0,length(rand_idx),1)
temp<-as.vector(biased_probability[1:3,1:3])
temp<-temp[temp!=0]
a<-sort(temp, index.return=TRUE)
temp_2<-a$x
temp_cum<-cumsum(temp_2)
idx[which(rand_idx<temp_cum[1])]=1
idx[which(temp_cum[1]<rand_idx & rand_idx<temp_cum[2])]=2
idx[which(temp_cum[2]<rand_idx & rand_idx<temp_cum[3])]=3
idx[which(temp_cum[3]<rand_idx)]=4

Sample<-matrix(0,n,2)
for(i in 1:4)
{
  switch(toString(a$ix[i]),
  '1'= Sample[idx==1,]<-matrix(t(rep(c(0,1),length(which(idx==1)))), length(which(idx==1)), byrow=TRUE),
  '2'= Sample[idx==2,]<-matrix(t(rep(c(1,0),length(which(idx==2)))), length(which(idx==2)), byrow=TRUE),
  '3'= Sample[idx==3,]<-matrix(t(rep(c(0,0),length(which(idx==3)))), length(which(idx==3)), byrow=TRUE),
  '4'= Sample[idx==4,]<-matrix(t(rep(c(1,1),length(which(idx==4)))), length(which(idx==4)), byrow=TRUE))}
#############################
#4. Define the weight matrix:
Weights_Matrix_1<-matrix(0, dim(Sample)[1],dim(Sample)[1])
Weights_Matrix_independent<-matrix(0, dim(Sample)[1],dim(Sample)[1])
for(i in seq(1,dim(Sample)[1],1))
{
  temp<-cbind(rep(Sample[i,1], dim(Sample)[1]), Sample[,2])  
  Weights_Matrix_1[i,]<-biased_probability[Sample[i,1]+1, temp[,2]+1]
  Weights_Matrix_independent[i,]<-ifelse(biased_probability[Sample[i,1]+1, temp[,2]+1]>0, 1,0)
}

#5. Sample valid permutations:
TargetSampleSize = 10
Permutations<-MCMC_Permutations(NULL,Weights_Matrix_1,TargetSampleSize,dim(Weights_Matrix_1)[1])
Permutations_independent<-MCMC_Permutations(NULL,Weights_Matrix_independent,TargetSampleSize,dim(Weights_Matrix_independent)[1])
#6. We are set:
temp<-matrix(0,TargetSampleSize,1)
temp_independent<-matrix(0,TargetSampleSize,1)
####################################
p_test<-c(1,0)#c(1,1)#c(2,2)#c(2,0)#c(1,0)

for(i in seq(1,TargetSampleSize,1))
{
  permuted_data<-cbind(Sample[,1], Sample[Permutations[,i],2])
  temp[i]<-length(which(permuted_data[,1]==p_test[1] & permuted_data[,2]==p_test[2] ))/dim(permuted_data)[1]
}
cat('Sampling with the induced distribution: ', toString(mean(temp)))

for(i in seq(1,TargetSampleSize,1))
{
  permuted_data<-cbind(Sample[,1], Sample[Permutations_independent[,i],2])
  temp_independent[i]<-length(which(permuted_data[,1]==p_test[1] & permuted_data[,2]==p_test[2] ))/dim(permuted_data)[1]
}
cat('Sampling with the P_W: ', toString(mean(temp_independent)))
###################################################################
W<-matrix(1,3,3)
W[3,]<-0
W[,3]<-0

P_w_x<-c(0.5,0.5,0)
P_w_y<-c(0.5,0.5,0)
#######
g_x<-c(0.5, 0.25, 0.25)

for(i in seq(1,num_of_repetitions,1))
{
  g_y<-P_w_y/(g_x*)
}