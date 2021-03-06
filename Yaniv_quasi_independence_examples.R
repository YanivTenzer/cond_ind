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
S<-sum(colSums(probability))-(sum(probability[1,3])+sum(probability[2,2]))
#Normalize the resulting distribution:
biased_probability<-probability/S
biased_probability[2,2]=0
biased_probability[1,3]=0
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
idx[which(temp_cum[3]<rand_idx & rand_idx<temp_cum[4])]=4
idx[which(temp_cum[4]<rand_idx & rand_idx<temp_cum[5])]=5
idx[which(temp_cum[5]<rand_idx & rand_idx<temp_cum[6])]=6
idx[which(temp_cum[6]<rand_idx)]=7
Sample<-matrix(0,n,2)
for(i in 1:7)
{
  switch(toString(a$ix[i]),
         '1'= Sample[idx==1,]<-matrix(t(rep(c(1,0),length(which(idx==1)))), length(which(idx==1)), byrow=TRUE),
         '2'= Sample[idx==2,]<-matrix(t(rep(c(2,0),length(which(idx==2)))), length(which(idx==2)), byrow=TRUE),
         '3'= Sample[idx==3,]<-matrix(t(rep(c(0,1),length(which(idx==3)))), length(which(idx==3)), byrow=TRUE),
         '4'= Sample[idx==4,]<-matrix(t(rep(c(2,1),length(which(idx==4)))), length(which(idx==4)), byrow=TRUE),
         '5'= Sample[idx==5,]<-matrix(t(rep(c(1,2),length(which(idx==5)))), length(which(idx==5)), byrow=TRUE),
         '6'= Sample[idx==6,]<-matrix(t(rep(c(0,0),length(which(idx==6)))), length(which(idx==6)), byrow=TRUE),
         '7'= Sample[idx==7,]<-matrix(t(rep(c(2,2),length(which(idx==7)))), length(which(idx==7)), byrow=TRUE))}
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
p_test<-c(2,0)#c(1,1)#c(2,2)#c(2,0)#c(1,0)

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
#####################################################################
cat("\014")
rm(list=ls())
path = 'D:/cond_ind_local'
setwd(path)
source('MCMC_Permutations.R')  
########################
#This is a toy example, simulating bivariate vector (x,y) \in {0,1,2}x{0,1,2}):
#1. Define P(X,Y):
probability<-matrix(0,2,2)
probability[1,1]=1/2
probability[2,2]=1/4
probability[1,2]=1/8
probability[2,1]=1/8
#2. W*P(x,y):
biased_probability<-matrix(0,2,2)
S<-sum(colSums(probability))-(sum(probability[2,2]))
#Normalize the resulting distribution:
biased_probability<-probability/S
biased_probability[2,2]=0
#3. Sample from  W*P(x,y):
n=10000
rand_idx<-runif(n)
idx<-matrix(0,length(rand_idx),1)
temp<-as.vector(biased_probability[1:2,1:2])
temp<-temp[temp!=0]
a<-sort(temp, index.return=TRUE)
temp_2<-a$x
temp_cum<-cumsum(temp_2)
idx[which(rand_idx<temp_cum[1])]=1
idx[which(temp_cum[1]<rand_idx & rand_idx<temp_cum[2])]=2
idx[which(temp_cum[2]<rand_idx)]=3
Sample<-matrix(0,n,2)
for(i in 1:3)
{
  switch(toString(a$ix[i]),
         '1'= Sample[idx==1,]<-matrix(t(rep(c(2,1),length(which(idx==1)))), length(which(idx==1)), byrow=TRUE),
         '2'= Sample[idx==2,]<-matrix(t(rep(c(1,2),length(which(idx==2)))), length(which(idx==2)), byrow=TRUE),
         '3'= Sample[idx==3,]<-matrix(t(rep(c(1,1),length(which(idx==3)))), length(which(idx==3)), byrow=TRUE))}
#############################
#4. Define the weight matrix:
Weights_Matrix_1<-matrix(0, dim(Sample)[1],dim(Sample)[1])
Weights_Matrix_independent<-matrix(0, dim(Sample)[1],dim(Sample)[1])
for(i in seq(1,dim(Sample)[1],1))
{
  temp<-cbind(rep(Sample[i,1], dim(Sample)[1]), Sample[,2])  
  Weights_Matrix_1[i,]<-biased_probability[Sample[i,1], temp[,2]]
  Weights_Matrix_independent[i,]<-ifelse(biased_probability[Sample[i,1], temp[,2]]>0, 1,0)
}

#5. Sample valid permutations:
TargetSampleSize = 10
Permutations<-MCMC_Permutations(NULL,Weights_Matrix_1,TargetSampleSize,dim(Weights_Matrix_1)[1])
Permutations_independent<-MCMC_Permutations(NULL,Weights_Matrix_independent,TargetSampleSize,dim(Weights_Matrix_independent)[1])
#6. We are set:
temp<-matrix(0,TargetSampleSize,1)
temp_independent<-matrix(0,TargetSampleSize,1)
####################################
p_test<-c(1,2)#c(1,1)#c(2,2)#c(2,0)#c(1,0)

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
#This is a toy example, simulating bivariate vector (x,y) \in {0,1,2}x{0,1,2}):
#1. Define P(X,Y):
probability<-matrix(0,2,2)
probability[1,1]=2/5
probability[2,2]=1/5#2/5
probability[1,2]=3/10
probability[2,1]=1/10
#2. W*P(x,y):
biased_probability<-matrix(0,2,2)
S<-sum(colSums(probability))-(sum(probability[2,2]))
#Normalize the resulting distribution:
biased_probability<-probability/S
biased_probability[2,2]=0
#3. Sample from  W*P(x,y):
n=20000
rand_idx<-runif(n)
idx<-matrix(0,length(rand_idx),1)
temp<-as.vector(biased_probability[1:2,1:2])
temp<-temp[temp!=0]
a<-sort(temp, index.return=TRUE)
temp_2<-a$x
temp_cum<-cumsum(temp_2)
idx[which(rand_idx<temp_cum[1])]=1
idx[which(temp_cum[1]<rand_idx & rand_idx<temp_cum[2])]=2
idx[which(temp_cum[2]<rand_idx)]=3
Sample<-matrix(0,n,2)
for(i in 1:3)
{
  switch(toString(a$ix[i]),
         '1'= Sample[idx==1,]<-matrix(t(rep(c(2,1),length(which(idx==1)))), length(which(idx==1)), byrow=TRUE),
         '2'= Sample[idx==2,]<-matrix(t(rep(c(1,2),length(which(idx==2)))), length(which(idx==2)), byrow=TRUE),
         '3'= Sample[idx==3,]<-matrix(t(rep(c(1,1),length(which(idx==3)))), length(which(idx==3)), byrow=TRUE))}
#############################
#4. Define the weight matrix:
Weights_Matrix_1<-matrix(0, dim(Sample)[1],dim(Sample)[1])
Weights_Matrix_independent<-matrix(0, dim(Sample)[1],dim(Sample)[1])
for(i in seq(1,dim(Sample)[1],1))
{
  temp<-cbind(rep(Sample[i,1], dim(Sample)[1]), Sample[,2])  
  Weights_Matrix_1[i,]<-biased_probability[Sample[i,1], temp[,2]]
  Weights_Matrix_independent[i,]<-ifelse(biased_probability[Sample[i,1], temp[,2]]>0, 1,0)
}

#5. Sample valid permutations:
TargetSampleSize = 10
Permutations<-MCMC_Permutations(NULL,Weights_Matrix_1,TargetSampleSize,dim(Weights_Matrix_1)[1])
Permutations_independent<-MCMC_Permutations(NULL,Weights_Matrix_independent,TargetSampleSize,dim(Weights_Matrix_independent)[1])
#6. We are set:
temp<-matrix(0,TargetSampleSize,1)
temp_independent<-matrix(0,TargetSampleSize,1)
####################################
p_test<-c(2,1)#c(1,1)#c(2,2)#c(2,0)#c(1,0)

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
#This is a toy example, simulating bivariate vector (x,y) \in {0,1,2}x{0,1,2}):
#1. Define P(X,Y):
probability<-matrix(0,2,2)
probability[1,1]=2/5
probability[2,2]=2/5#2/5
probability[1,2]=1/10
probability[2,1]=1/10
#2. W*P(x,y):
biased_probability<-matrix(0,2,2)
S<-sum(colSums(probability))-(sum(probability[2,2]+probability[1,1]))
#Normalize the resulting distribution:
biased_probability<-probability/S
biased_probability[2,2]=0
biased_probability[1,1]=0
#3. Sample from  W*P(x,y):
n=20000
rand_idx<-runif(n)
idx<-matrix(0,length(rand_idx),1)
temp<-as.vector(biased_probability[1:2,1:2])
temp<-temp[temp!=0]
a<-sort(temp, index.return=TRUE)
temp_2<-a$x
temp_cum<-cumsum(temp_2)
idx[which(rand_idx<temp_cum[1])]=1
idx[which(temp_cum[1]<rand_idx & rand_idx<temp_cum[2])]=2
idx[which(temp_cum[2]<rand_idx)]=3
Sample<-matrix(0,n,2)
for(i in 1:2)
{
  switch(toString(a$ix[i]),
         '1'= Sample[idx==1,]<-matrix(t(rep(c(1,2),length(which(idx==1)))), length(which(idx==1)), byrow=TRUE),
         '2'= Sample[idx==2,]<-matrix(t(rep(c(2,1),length(which(idx==2)))), length(which(idx==2)), byrow=TRUE))}
#############################
#4. Define the weight matrix:
Weights_Matrix_1<-matrix(0, dim(Sample)[1],dim(Sample)[1])
Weights_Matrix_independent<-matrix(0, dim(Sample)[1],dim(Sample)[1])
for(i in seq(1,dim(Sample)[1],1))
{
  temp<-cbind(rep(Sample[i,1], dim(Sample)[1]), Sample[,2])  
  Weights_Matrix_1[i,]<-biased_probability[Sample[i,1], temp[,2]]
  Weights_Matrix_independent[i,]<-ifelse(biased_probability[Sample[i,1], temp[,2]]>0, 1,0)
}

#5. Sample valid permutations:
TargetSampleSize = 10
Permutations<-MCMC_Permutations(NULL,Weights_Matrix_1,TargetSampleSize,dim(Weights_Matrix_1)[1])
Permutations_independent<-MCMC_Permutations(NULL,Weights_Matrix_independent,TargetSampleSize,dim(Weights_Matrix_independent)[1])
#6. We are set:
temp<-matrix(0,TargetSampleSize,1)
temp_independent<-matrix(0,TargetSampleSize,1)
####################################
p_test<-c(2,1)#c(1,1)#c(2,2)#c(2,0)#c(1,0)

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
