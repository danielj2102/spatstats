# 170009243
# Daniel McGrogan

library(mvtnorm)
library(spatstat)


#simulations:
#===================
my_rThomas <- function(xdim,ydim,kappa,mu,scale){
  {
   #Extended simulation windows parameters
    extension <- 4*scale
    xminnew <- xdim[1]-extension
    xmaxnew <- xdim[2]+extension
    yminnew <- ydim[1]-extension
    ymaxnew <- ydim[2]+extension
    
    #new dimensions
    win_dimnewx <- xmaxnew-xminnew
    win_dimnewy <- ymaxnew-yminnew
    win_areanew <- win_dimnewx*win_dimnewy
    
    #simulate for parents
    no_ppoints <- rpois(1, win_areanew*kappa)
    
    #coordinates for parent points
    xparent <- xminnew+win_dimnewx*runif(no_ppoints)
    yparent <- yminnew+win_dimnewy*runif(no_ppoints)
    
    #Poisson point process for daughter points
    no_dpoints <- rpois(no_ppoints,mu)
    total_p <- sum(no_dpoints)
    
    #Generate the locations 
    #by simulating independent normal variables
    xloc <- rnorm(total_p,0,scale)
    yloc <- rnorm(total_p,0,scale)
    
    #centre of clusters
    xcent <- rep(xparent,no_dpoints)
    ycent <- rep(yparent,no_dpoints)
    
    #translate 
    x <- xcent+xloc
    y <- ycent+yloc
    
    #thin points if they are outside of the simulation window
    inside <- ((x>=xdim[1])&(x<=xdim[2])&(y>=ydim[1])&(y<=ydim[2]))
    
    #new points:
    x <- x[inside]
    y <- y[inside]
  }
  return(ppp(x,y))
}

set.seed(1234) #For reproducibility

#2 simulations:
sim1 <- my_rThomas(c(0,1),c(0,1),10,50,0.01)
sim2 <- my_rThomas(c(0,1),c(0,1),10,50,0.02)

#plot the 2 simulations:
plot(sim1,type= 'p',xlab='x',ylab='y')
plot(sim2,type='p',xlab='x',ylab='y')

#simulation one model fitting:

#Thomas model fitting: Simulation 1

sim1k <- Kest(sim1, correction = "translate") #k-function for simulation 1

thomas1 <- thomas.estK(sim1k, startpar = c(kappa=10,scale=0.01))
plot(thomas1)

thomas1$par # returns estimates of Kappa and sigma**2
estmu1 <- sim1$n/thomas1$par[1] #estimates mu

thomas.sim1 <- rThomas(kappa = thomas1$par[1], scale = sqrt(thomas1$par[2]), mu = estmu1)
plot(thomas.sim1, main = "sim 1 thomas")
plot(pcf(thomas.sim1)) #returns pcf of the thomas simulation with calculated parameters

#Matern model fitting: Simulation 1

mat1 <- matclust.estK(sim1k, startpar = c(kappa=10,scale=0.01))
plot(mat1)

mat1$par # returns estimates of Kappa and sigma
estmu2 <- sim1$n/mat1$par[1] #estimates mu

mat.sim1 <- rMatClust(kappa = mat1$par[1], scale = mat1$par[2], mu = estmu2)
plot(mat.sim1, main = "sim 1 matern")
plot(pcf(mat.sim1)) #returns pcf of the matern simulation with calculated parameter

#simulation two model fitting:

#Thomas model fitting: Simulation 2
sim2k <- Kest(sim2, correction = "translate") #k-function for simulation 1

thomas2 <- thomas.estK(sim2k, startpar = c(kappa=10,scale=0.02))
plot(thomas2)

thomas2$par # returns estimates of Kappa and sigma**2
estmu3 <- sim2$n/thomas2$par[1] #returns estimate of mu

thomas.sim2 <- rThomas(kappa = thomas2$par[1], scale = sqrt(thomas2$par[2]), mu= estmu3)
plot(thomas.sim2, mu = estmu3, main = "sim 2 thomas ")
plot(pcf(thomas.sim2)) #returns pcf of the thomas simulation with calculated parameter


#Matern model fitting: Simulation 1
mat2 <- matclust.estK(sim2k, startpar = c(kappa=10,scale=0.02))
plot(mat2)

mat2$par # returns estimates of Kappa and sigma
estmu4 <- sim2$n/mat2$par[1] #returns estimate of mu

mat.sim2 <- rMatClust(kappa = mat2$par[1], scale = mat2$par[2], mu = estmu4)
plot(mat.sim2, main = "sim 2 matern")
plot(pcf(mat.sim2)) #returns pcf of the matern simulation with calculated parameter


#fitting to data
#===================
attach(bei)
attach(bei.extra)

set.seed(1234) #For reproducibility

bw.dens<- bw.diggle(bei)
dens <- density.ppp(bei, bw.diggle) #initial estimate of intensity
plot(dens) #default kernal is gaussian

#alternative bandwidth selection method
bw.dens2 <- bw.scott(bei)
dens2 <- density.ppp(bei, bw.scott) #drastically different intensity
plot(dens2) 

#fitting model x&y coordinates, polynomial variants and interaction terms
xymodel <- ppm(bei~x+y)
summary(xymodel)

polyx2y <- ppm(bei~ polynom(x,2)+y)
#summary(polyx2y)

polyxy2 <- ppm(bei~ x+polynom(y,2))
#summary(polyxy2)

intmodel1 <- ppm(bei~ x+y+x*y)
#summary(intmodel1)

intmodel2 <- ppm(bei~ polynom(x,2)+polynom(y,2)+x*y)
#summary(intmodel2)

#compare AIC:
AIC(xymodel)
AIC(polyx2y)
AIC(polyxy2)
AIC(intmodel1)
AIC(intmodel2)

# Hence intmodel2 has the lowest AIC value and should be chosen
plot(intmodel2, se=F)

#gradient, altitude and interaction of these
grad <- bei.extra$grad
alt <- bei.extra$elev

gradaltmodel <- ppm(bei~grad+alt)
#summary(gradaltmodel)

intmodel3 <- ppm(bei~grad*alt)
#summary(intmodel3)

intmodel4 <- ppm(bei~ polynom(grad,2)+alt+polynom(grad*alt,2))
#summary(intmodel4)

intmodel5 <- ppm(bei~ grad+polynom(alt,2)+polynom(grad*alt,2))
#summary(intmodel5)

#compare AIC:
AIC(gradaltmodel)
AIC(intmodel3)
AIC(intmodel4)
AIC(intmodel5) #select intmodel5 as AIC lowest

plot(intmodel5, se=F)

intmodel5pcf <- pcfinhom(bei, intmodel5)
plot(intmodel5pcf) #pcf for inhom model
intmodel5k <- Kinhom(bei,intmodel5,correction=c("translate", "Ripley"),nlarge = 5000)
plot(intmodel5k)#k-function for inhom model

#preffered model structure: intmodel5
coxintmodel5 <- kppm(bei~ grad+polynom(alt,2)+polynom(grad*alt,2), clusters = "Thomas")
#summary(coxintmodel5)
plot(coxintmodel5)
par(ask=F)
par(ask=F) #bypasses "Hit <Return> to see next plot"
coxintmodel5pcf <- pcfinhom(bei, coxintmodel5)
plot(coxintmodel5pcf) #pcf for cox model
coxintmodel5k <- Kinhom(bei,coxintmodel5,correction=c("translate", "Ripley"),nlarge = 5000)
plot(coxintmodel5k) #k-function for cox model


set.seed(1234) #for reproducibility
coxsim <-simulate.kppm(coxintmodel5,nsim=2)
par(mfrow=c(2,1))
plot(coxsim) #2 simulations from the Cox process


