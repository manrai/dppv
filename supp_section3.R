library(poibin)
library(tidyverse)
library(ggthemes)
library(plotly)
library(reshape2)

naive <- function(pp,iter) {
  P = numeric(iter)
  for (j in 1:iter) {
    s = 0
    for (i in 1:length(pp)) {
      s <- s + rbinom(1,1,pp[i])
    }
    P[j] = s
  }
  return(P)
}
dftcf <- function(pp,iter) {
  return(rpoibin(iter,pp))
}
rna_p <- function(n,iter) {
  P = numeric(iter)
  for (i in 1:iter) {
    pp = runif(n)
    P <- mean(ppoibin(kk=0:n, pp=pp, method = "RNA",wts=NULL))
  }
  return(P)
}
dftcf_p <- function(n,iter) {
  P = numeric(iter)
  for (i in 1:iter) {
    pp = runif(n)
    P <- mean(ppoibin(kk=0:n, pp=pp, method = "DFT-CF",wts=NULL))
  }
  return(P)
}
discrete.inv.transform.sample <- function( p.vec ) {
  U  <- runif(1)
  if(U <= p.vec[1]){
    return(1)
  }
  for(state in 2:length(p.vec)) {
    if(sum(p.vec[1:(state-1)]) < U && U <= sum(p.vec[1:state]) ) {
      return(state)
    }
  }
}

##############################
###### Runtime Analysis ######
##############################

# Random variate generation
# Naive
system.time(naive(pp=runif(500),iter=100)) # 0.1 seconds
system.time(naive(pp=runif(5000),iter=100)) # 1 second
system.time(naive(pp=runif(50000),iter=100)) # 10 seconds

system.time(naive(pp=runif(5000),iter=10)) # 0.1 seconds
system.time(naive(pp=runif(5000),iter=100)) # 1 second
system.time(naive(pp=runif(5000),iter=1000)) # 10 seconds

# DFT-CF
system.time(dftcf(pp=runif(500),iter=100)) # 0.005 seconds
system.time(dftcf(pp=runif(5000),iter=100)) # 0.352 second
system.time(dftcf(pp=runif(50000),iter=100)) # 34 seconds

system.time(dftcf(pp=runif(5000),iter=10)) # 0.3 seconds
system.time(dftcf(pp=runif(5000),iter=100)) # 0.3 seconds
system.time(dftcf(pp=runif(5000),iter=1000)) # 0.3 seconds

# CDF generation
#DFT-CF
system.time(dftcf_p(n=500,iter=10))
system.time(dftcf_p(n=5000,iter=10))
system.time(dftcf_p(n=50000,iter=10))

# RNA
system.time(rna_p(n=500,iter=100))
system.time(rna_p(n=5000,iter=100))
system.time(rna_p(n=50000,iter=100))

##############################
###### Approx. Analysis ######
##############################

iter = 1000
n = 100
num.variates = 5000
m_exp = numeric(iter)
m_obs = numeric(iter)
v_exp = numeric(iter)
v_obs = numeric(iter)

for (i in 1:iter) {
  a1 = runif(1)
  b1 = runif(1)
  pp = rbeta(n,a1,b1)
  a2 = runif(1)
  b2 = runif(1)
  qq = rbeta(n,a2,b2)
  
  p = rpoibin(num.variates,pp=pp)
  q = rpoibin(num.variates,pp=qq)
  dppv = p/(p+q)
  m1 = mean(dppv) # "exact"
  v1 = var(dppv) # "exact"
  m2 = mean_est = sum(pp) / (sum(pp)+sum(qq)) # delta method
  v2 = var_est = (sum(qq)^2 * var(p) + sum(pp)^2 * var(q)) / (sum(pp) + sum(qq))^4 # delta method
  
  m_exp[i] = m1
  m_obs[i] = m2
  v_exp[i] = v1
  v_obs[i] = v2
  print(i)
}

# pdf(file="plots/mean_var_check.pdf",width = 8,height=5)
# par(mfrow=c(1,2))
# plot(m_exp,m_obs,xlim = c(0,1),ylim=c(0,1))
# abline(a=0,b=1,col='red')
# plot(v_exp,v_obs,xlim=c(0,max(c(v_exp,v_obs),na.rm = T)),ylim=c(0,max(c(v_exp,v_obs),na.rm = T)))
# abline(a=0,b=1,col='red')
# dev.off()

m_diffs <- abs(m_exp-m_obs)/m_exp
v_diffs <- abs(v_exp-v_obs)/v_exp
100*mean(m_diffs[is.finite(m_diffs)])
100*mean(v_diffs[is.finite(v_diffs)])

##############################
### Homogeneity of the p_i ###
##############################

# beta distribution plot
alpha = seq(0,3,by=0.01)
beta = seq(0,3,by=0.01)
a.b <- as.tibble(merge(alpha,beta))
colnames(a.b) <- c('a','b')
a <- a.b$a
b <- a.b$b
a.b$var <- (a*b)/( ((a+b)^2) * (a+b+1) )

a.b.wide <- data.matrix(acast(a.b,a~b,value.var='var'))

axx <- list(
  title = "alpha"
)

axy <- list(
  title = "beta"
)

axz <- list(
  title = "variance"
)

p <- plot_ly(x = ~a, y = ~b, z = ~a.b$var, type = "mesh3d",
             colorscale = 'Viridis',alpha = 0.5) %>%
  layout(scene = list(xaxis=axx,yaxis=axy,zaxis = axz))
p

# contour version
Variance = a.b$var
f2 <- list(
  family = "Arial, sans-serif",
  size = 18,
  color = "black"
)
aaa <- list(
  title = "",
  showticklabels = TRUE,
  tickangle = -45,
  tickfont = f2,
  exponentformat = "E"
)
l <- list(
  title="Variance"
)
p <- plot_ly(x = ~a, y = ~b, z = ~Variance, type = "contour") %>% layout(xaxis=aaa,yaxis=aaa,legend=l)
p

# heatmap version
p <- plot_ly(x = ~a, y = ~b, z = ~a.b$var, type = "heatmap")
p

# variance vs. beta value plot
# consider p_i ~ beta(alpha,beta=alpha) --> variance is 1/(4*(2*beta+1))
iter = 100
beta.vec = seq(0,3,(3/iter))
var.pp = 1/(4*(2*beta.vec+1))
num.variates = 5000
n = 500
var.vec = numeric(iter)
qq = rep(0.5,500)

for (i in seq_along(beta.vec)) {
  beta = beta.vec[i]
  pp = rbeta(n,shape1=beta,shape2=beta)
  p = rpoibin(m = num.variates, pp)
  q = rpoibin(m = num.variates, qq)
  dppv = p/(p+q)
  var.vec[i] <- var(dppv)
}

df <- tibble(beta.vec,var.vec,var.pp)
p <- ggplot(df,aes(beta.vec)) +
  geom_point(aes(y = var.vec), colour = 'red') +
  geom_smooth(aes(y = var.vec)) +
  theme_tufte()
p
ggsave(filename = "plots/var_vs_beta.pdf",plot=p)

##############################
########## dPPV & n ##########
##############################

iter = 10
mean.vec <- numeric(iter)
var.vec <- numeric(iter)
num.variates = 5000
n = 500
pp = rbeta(500,2,2)
qq = rep(0.5,500)

for (i in 1:iter) {
  p = rpoibin(m = num.variates, rep(pp,i))
  q = rpoibin(m = num.variates, rep(qq,i))
  dppv = p/(p+q)
  mean.vec[i] <- mean(dppv,na.rm=TRUE)
  var.vec[i] <- var(dppv,na.rm=TRUE)
}

df <- tibble(k=1:iter,m=mean.vec,v=1000*var.vec)
p <- ggplot(df,aes(k)) +
  geom_point(aes(y=m),colour='red') +
  geom_point(aes(y=v),colour='blue') +
  theme_tufte()
p
ggsave(filename="plots/dppv_vs_n.pdf",plot=p)





