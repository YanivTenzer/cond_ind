library('mvtnorm')
rho=0.6
sigma<-matrix(c(1, rho, rho , 1),2,2)
mu=c(0,0)
n=10^4
data<-rmvnorm(n, mu, sigma)

truncated_data<-data[data[,2]>data[,1],]
X_w_cdf<-ecdf(truncated_data[,1])
Y_w_cdf<-ecdf(truncated_data[,2])

u<-runif(n=10000)
v<-runif(n=10000)

X_w_Y_w<-cbind(quantile(truncated_data[,1], u),quantile(truncated_data[,2], v))
WX_w_Y_w<-X_w_Y_w[X_w_Y_w[,2]>X_w_Y_w[,1],]
###############################
library('kdecopula')
library('copula')

gumbel.cop <- gumbelCopula(3)
x <- rCopula(400, gumbel.cop)
true_density<-dCopula(x,gumbel.cop)

dens.est <- kdecop(x)
rho=cor(x,method="spearman")[1,2]
n=dim(x)[1]
BW=max(1, round(n^(1/3) * exp(abs(rho)^(1/n)) * (abs(rho) + 0.1)))
kdecopObj<-kdecop(x, bw = 65, mult = 1, method = "bern", knots = 30)
density_est<-dkdecop(x,kdecopObj , stable = TRUE)

temp<-cbind(true_density,density_est)