

#######################################################################
#   simulate  cyton model
#######################################################################
#calculations based on explicit analytic formulas
#total number of dividing cell; an explicit expression for r_n^div (defined implicity
#in eqn 13 of Miao et al) was derived and then inserted into eqn 15 of Miao et al to get
#the total number of cells.
sum.rn.div2 <- function(i,a,n,t,theta){
  j <- 0:i
  return ((-1)^(n+i+1)/a^(n-i)*(1/theta^(i+1) - exp(-theta*t)*sum(t^(i-j)/factorial(i-j)/theta^(j+1))))
}
sum.rn.div <- function(n,t,par){
  theta0 <- 1/par['mu0p'] + 1/par['mu0d']
  if(n == 0){
    return (par['f0']*N/par['mu0p']/theta0*(1 - exp(-theta0*t)))
  }
  else{
    theta <- 1/par['mup'] + 1/par['mud']
    a <- theta0 - theta
    temp <- par['f0']/par['mu0p']*(2*par['f1']/par['mup'])^n*N
    temp <- temp*(((-1/a)^n*(1 - exp(-theta0*t))/theta0) +
                    sum(sapply(0:(n-1),sum.rn.div2,a,n,t,theta)))
    return (temp)
  }
}
#calculations based on explicit analytic formulas
#total number of dieing cells; an explicit expression for r_n^die (defined implicity
#in eqn 14 of Miao et al) was derived and then inserted into eqn 15 of Miao et al to get
#the total number of cells.
sum.rn.die2 <- function(i,a,n,t,theta){
  j <- 0:(i+1)
  return ((-1)^(n+i)/a^(n-i-1)*(1/theta^(i+2) - exp(-theta*t)*sum(t^(i-j+1)/factorial(i-j+1)/theta^(j+1))))
}
sum.rn.die <- function(n,t,par){
  theta0 <- 1/par['mu0p'] + 1/par['mu0d']
  if(n == 0){
    return (N/par['mu0d']*((1-par['f0'])*par['mu0d']*(1 - exp(-t/par['mu0d'])) +
                             par['f0']/theta0*(1 - exp(-theta0*t))))
  }
  else{
    theta <- 1/par['mup'] + 1/par['mud']
    a <- theta0 - theta; b <- 1/par['mu0d'] - theta
    temp <- 2*N*par['f0']/par['mu0p']/par['mud']*(2*par['f1']/par['mup'])^(n-1)
    temp2 <- (-1)^(n-1)/a^n*((1 - exp(-theta*t))/theta + (exp(-theta0*t) - 1)/theta0) +
      sum(sapply(0:(n-2),sum.rn.die2,a,n,t,theta))*(n > 1)
    return (temp*temp2)
  }
}
temp <- rep(0:ngen,each=length(t_pts))
sim.data <- mapply(sum.cells,split(temp,1:length(temp)),split(rep(t_pts,ngen+1),1:length(temp)),rep(list(parms),length(temp)))
sim.data <- matrix(sim.data,ncol=ngen+1,byrow=FALSE)
print(sim.data)
(sum((sim.data - real.data)^2))

#this function implements eqn (15) of the Miao et al by calling the other functions
#defined above
sum.cells <- function(n,t,par){
  return(ifelse(n == 0, N - sum.rn.div(n,t,par) - sum.rn.die(n,t,par),
                2*sum.rn.div(n-1,t,par) - sum.rn.div(n,t,par) - sum.rn.die(n,t,par)))
}


##################################################################
# Optimization of cyton model
##################################################################
ngen <- 20
t_pts <- seq(0,120,12)          # initialising default parameters
N <- 50000           
###############   parameter for cyton model #######################

# read the simulated data generated from agent based model
df <- read.table("./datasets/agent_based_sim_output_3.txt",header = TRUE)

sum.cols <- function(idx, data) return (rowSums(data[,idx]))
# reshape data to have total number of cells in each division class in a single column
real.data <- matrix( rep(0, (length(t_pts)*(ngen+1))), nrow = length(t_pts))
real.data[1,1] <- N[1]
real.data[-1,] <- sapply(split(2:43,rep(1:21,each=2)),sum.cols,df)

# this function takes as inputs a parameter vector (called parms)
#and a "real" dataset (called real.data). It generates a simulated data
#(called sim.data) corresponding to input parameters and the calculates
#the sum of squared difference (SSD) between real.data and sim.data
#generates 
fit.func <- function(parms, real.data){
  temp <- rep(0:ngen,each=length(t_pts))
  sim.data <- mapply(sum.cells,split(temp,1:length(temp)),split(rep(t_pts,ngen+1),1:length(temp)),rep(list(parms),length(temp)))
  sim.data <- matrix(sim.data,ncol=ngen+1,byrow=FALSE)
  return(sum((sim.data - real.data)^2))
}
# optim function determines the best parameter canditate
lb <- rep(1e-3,6) #lower-bound of paramter values
ub <- c(rep(1,4),rep(4,100)) #upper-bound of parameter values
parms <- runif(length(lb),min=lb,max=ub) #new random parameter vector
names(parms) <- c('f0', 'f1', 'mu0d', 'mud', 'mu0p','mup') #names of the parameters
opt.parms <- optim(parms, fit.func, gr=NULL, real.data)
print(opt.parms)
