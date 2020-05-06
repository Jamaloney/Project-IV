#Masterfile of code
#section 1: Introduction
#There is no R code presented in this chapter

#Section 2: 
#Spectral densities
par(mfrow=c(1,2))

curve(dnorm(x, mean=0, sd=0.7), from = -6, to = 6, col = "#F8766D", xlab = "x", ylab = "Probability Density", lwd = 2)
curve(dnorm(x, mean=0, sd=1), from = -6, to = 6, col = "#00BFC4", add = TRUE, lwd = 2)
curve(dnorm(x, mean=0, sd=3), from = -6, to = 6, col = "#7CAE00", add = TRUE, lwd = 2)
curve(dnorm(x, mean=0, sd=5), from = -6, to = 6, col = "#C77CFF", add = TRUE, lwd = 2)
legend("topleft", legend = paste0("σ = ", c(0.7,1,3,5)),
       col = c("#F8766D","#00BFC4","#7CAE00","#C77CFF"),
       lty = 1, lwd = 2)

curve(dt(x, 30), from = -6, to = 6, col = "#F8766D", 
      xlab = "x", ylab = "Probability Density", lwd = 2)
curve(dt(x, 10), from = -6, to = 6, col = "#00BFC4", add = TRUE, lwd = 2)
curve(dt(x, 3), from = -6, to = 6, col = "#7CAE00", add = TRUE, lwd = 2)
curve(dt(x, 1), from = -6, to = 6, col = "#C77CFF", add = TRUE, lwd = 2)
legend("topleft", legend = paste0("m = ", c(1, 3, 10, 30)),
       col = c("#F8766D","#00BFC4","#7CAE00","#C77CFF"),
       lty = 1, lwd = 2)

#Regression model vs 1D emulations
#1D emulation function
#Code regarding 1D emulation and implausibility is based on code from a Durham emulation workshop

emulator_1d <- function(x_j,x,D,beta_0=0,sig2=1,theta=15){
  
  nx_j <- length(x_j)			
  nx <- length(x)					
  x_both <- c(x,x_j)			
  n <- nx + nx_j					
  
  E_fx <- rep(beta_0,nx)			
  E_D <- rep(beta_0,nx_j)			
  
  Sigma <- sig2 * exp( - as.matrix(dist(x_both))^2/theta^2 )   
  
  Var_fx 	 <- Sigma[1:nx,1:nx]
  Var_D 	 <- Sigma[(nx+1):n,(nx+1):n]
  Cov_fx_D <- Sigma[1:nx,(nx+1):n]							
  
  ED_fx <- E_fx + Cov_fx_D %*% solve(Var_D) %*% (D - E_D)							
  VarD_fx <- Var_fx - Cov_fx_D %*% solve(Var_D) %*% t(Cov_fx_D)					
  sdD_fx <- sqrt(diag(VarD_fx))
  return(list(Expect=ED_fx,StanDev=sdD_fx))
}

#Function to plot 1D emulation/implausibility
#green #blue #red #orange "#81d381","#6b9aea","#ff6961, #ffb347"
plot_emulator_1d_pastelcolours <- function(ED_fx,sdD_fx,z=-0.5,sigma_e=0.02,sigma_epsilon=0,yloc=-1.35,ylimit,plottype="emul_imp"){
  
  if(plottype=="emul"){
    plot(x,ED_fx,ty="l",col="#6b9aea",ylim=ylimit,lwd=2,xlab=expression(I[i:j]),ylab="Emulator Output")
    lines(x,ED_fx+2*sdD_fx,col="#ff6961",lwd=2)
    lines(x,ED_fx-2*sdD_fx,col="#ff6961",lwd=2)
    points(x_j,D,pch=16,cex=1.5)				
  }
  
  if(plottype=="imp"){
    Imp <- sqrt(  (ED_fx - z)^2 / ( sdD_fx^2 + sigma_epsilon^2 + sigma_e^2)  )
    plot(x,Imp,ty="l",ylim=c(0,30),xlab="Input Parameter x",ylab="Implausibility",lwd=2)
    abline(h=3,col="#6b9aea",lwd=2)
    
    points(x[Imp>5],rep(yloc,sum(Imp>5)),pch=16,col="#ff6961",cex=1.5)				
    points(x[(3<Imp) & (Imp<5)],rep(yloc,sum((3<Imp) & (Imp<5))),pch=16,col="#ffb347",cex=1.5)
    points(x[Imp<3],rep(yloc,sum(Imp<3)),pch=16,col="#81d381",cex=1.5)
  }
  
  if(plottype=="emul_imp"){
    
    plot(x,ED_fx,ty="l",col="grey40",ylim=ylimit,lwd=2,xlab="Input Parameter x",ylab="Emulator Output")
    polygon(c(min(x)-1,min(x)-1,max(x)+1,max(x)+1),c(z-3*sigma_e,z+3*sigma_e,z+3*sigma_e,z-3*sigma_e),col="#ffe6e5",border = NA)
    lines(x,ED_fx,col="grey40",lwd=2)
    lines(x,ED_fx+2*sdD_fx,col="grey60",lwd=2)
    
    lines(x,ED_fx-2*sdD_fx,col="gray60",lwd=2)
    points(x_j,D,pch=16,cex=1.5)				
    
    Imp <- sqrt(  (ED_fx - z)^2 / ( sdD_fx^2 + sigma_epsilon^2 + sigma_e^2)  )
    
    
    abline(h=c(z,z+3*sigma_e,z-3*sigma_e),lty=c(1,2,2),lwd=1.4)
    
    
    points(x[Imp>5],rep(yloc,sum(Imp>5)),pch=16,col="#ff6961",cex=1.5)				
    points(x[(3<Imp) & (Imp<5)],rep(yloc,sum((3<Imp) & (Imp<5))),pch=16,col="#ffb347",cex=1.5)
    points(x[Imp<3],rep(yloc,sum(Imp<3)),pch=16,col="#81d381",cex=1.5)
    print((x[Imp<3]))
  }
}

#Running an example
#Take data points x_j
x_j <- c(1,2,3,4)	
nx <- 200										#emulate at 200 points			
x <- seq(0.7,04.3,len=nx)					
fx_j <- 0.3*exp(((x_j)/2)) -1.5				
D <- fx_j			#D referred to as Z in presentation
plot(D,pch=16,cex=1.5,ylim=c(-1.3,1.3),lwd=2,xlab="Input Parameter x",ylab="Output f(x)")
out <- emulator_1d(x_j=x_j,x=x,D=D,beta_0=0,sig2=1,theta=1.5)

plot_emulator_1d_pastelcolours(ED_fx=out$Expect,sdD_fx=out$StanDev,z=-0.5,sigma_e=0.02,sigma_epsilon=0,yloc=-1.35,ylimit=c(-1.3,1.3),plottype="emul")
legend("topleft", legend=c("Adjusted expectation", "3SD credible interval" ),
       col=c("#6b9aea","#ff6961"), lty=c(1,1),lwd=c(2,2) ,cex=0.8)

#compare to a regression plot
plot(D,pch=16,cex=1.5,ylim=c(-1.3,1.3),lwd=2,xlab="Input Parameter x",ylab="Output h(x)")
abline(lm(fx_j~ x_j), lty = 2,lwd=4, col = "#81d381")
lm(fx_j~ x_j)
curve(0.3*exp(((x)/2)) -1.5	, add = TRUE, col = "black")
conf_interval <- predict(lm(fx_j~ x_j),  interval="confidence",level = 0.95)
lines(conf_interval[,2], col=2, lty=2, lwd=2.5)
lines(conf_interval[,3], col=2, lty=2,lwd=2.5)
legend("topleft", legend=c("True function", "Regression line", "3SD confidence interval"),
       col=c("black", "#81d381",2), lty=c(1,2,2),lwd=c(1,2,1.5), cex=0.8)


#Illustrations of basics of LHS
library(tgp)
lhsexmp <- lhs(15, rbind(c(0,5), c(-2.5,3.3 )))
plot(lhsexmp,xlab=expression(x[1]),ylab=expression(x[2]),pch=16,cex=0.7,col='blue')
abline(v=seq(0,5,length=16), lty=2, col=3)
abline(h=seq(-2.5,3.3,length=16), lty=2, col=3)



#Pi example 
circle_plotter = function(n){
  x=runif(n)
  y=runif(n)
  z=sqrt(x^2+y^2)
  length(which(z<=1))*4/length(z)
  curve(sqrt(1-x^2), from=0, to=1, , xlab="x", ylab="y",lwd=3)
  points(x[which(z<=1)],y[which(z<=1)],xlab="",ylab="",main=" ",col="blue",pch=19)
  points(x[which(z>1)],y[which(z>1)],col='coral',pch=19)
  
}

circle_plotter(100)
circle_plotter(2000)


LHS_method = function(long){
  c = rep(0,long)
  numberIn = 0
  x = randomLHS(long,2)
  for(i in 1:long){
    if(sqrt(x[i,1]^2 + x[i,2]^2) <= 1){
      numberIn = numberIn + 1
    }
    prop = numberIn / i
    piHat = prop *4
    c[i] = piHat
  }
  return(c)
}
aug_LHS_method = function(long){
  c = rep(0,long)
  numberIn = 0
  x = randomLHS(1,2)
  for(j in 1:long){
    x= augmentLHS(x,1)
  }
  for(i in 1:long){
    
    if(sqrt(x[i,1]^2 + x[i,2]^2) <= 1){
      numberIn = numberIn + 1
    }
    prop = numberIn / i
    piHat = prop *4
    c[i] = piHat
  }
  return(c)
}


#plot with augmented LHS and updated LHS
size1 = 600
res1 = simulation(size1)
res2 = aug_LHS_method(size1)
ini1 = 1

#Plot with PLHS, MC and LHS, 600 samples
plot(LHS_method(size1),type = "l",lty=1, xlab = "Number of samples", ylab = "Estimate of π")
lines(rep(pi, size1)[ini:size], col="red", lty = 4)
lines(res1[ini:size1],col=' dark grey')
lines(res2[ini:size1],col='#00BFC4')
legend("topright", legend = paste0(c("π","MC","LHS","PLHS")),
       col = c("red", "dark grey ", "  black","#00BFC4"),
       lty = c(4,1,1,1), lwd = 2)
#Plot with MC and LHS
size = 7000
res = simulation(size)
ini = 1
plot(simulation_LHS(size),type = "l",lty=1, xlab = "Number of samples", ylab = "Estimate of π")
lines(rep(pi, size)[ini:size], col="red", lty = 4)
lines(res[ini:size],col=' dark grey')
legend("topright", legend = paste0(c("π","MC","LHS")),
       col = c("red", "dark grey ", "  black"),
       lty = c(4,1,1), lwd = 2)


#2D emulation example
xrange <- rbind(c(-4,1),c(-4,1))				
nx_j <- 15
x_j <- lhs(n=nx_j,rect=xrange)				
plot(x_j,pch=16, xlab = expression(x[1]), ylab= expression(x[2]))							

#grid of points in 2d to evaluate the emulator at
x1seq <- seq(-4,1,len=15)
x2seq <- seq(-4,1,len=15)
x <- as.matrix(expand.grid(x1seq,x2seq))	

#run the model
fx_j <- 2*sin(-0.5*x_j[,1])+cos(1.5*x_j[,2])				
D <- fx_j								
#use the emulator to predict at all the points in the grid represented by x
emulator_2d <- function(x_j,x,D,beta_0=0,sig2=1,theta=c(1,1),delta=10^-5){
  
  
  nx_j <- nrow(x_j)				# number of points in x_j 
  nx <- nrow(x)					# number of points we want to evaluate the emulator for (apply the BL update for) 
  x_both <- rbind(x,x_j)			# join both x and x_j together (just to make covariance calculation easier)
  n <- nx + nx_j					# total number of old and new points
  #Prior expectation of f(x) and the vector D = f(x_j) 
  E_fx <- rep(beta_0,nx)			# needed for BL update
  E_D  <- rep(beta_0,nx_j)		# needed for BL update
  
  #Prior variances and covariances of f(x) and D = f(x_j)
  Sigma <- sig2 *( (1-delta)* exp( - as.matrix(dist( t( t(x_both)/theta) )^2 ) ) + diag(delta,n))  # bit that makes f(x) a smooth function!
    Var_fx 	 <- Sigma[1:nx,1:nx]
  Var_D 	 <- Sigma[(nx+1):n,(nx+1):n]
  Cov_fx_D <- Sigma[1:nx,(nx+1):n]
  
  #Bayes Linear Update equations.
  ED_fx <- E_fx + Cov_fx_D %*% solve(Var_D) %*% (D - E_D)						
  VarD_fx <- Var_fx - Cov_fx_D %*% solve(Var_D) %*% t(Cov_fx_D)			
  sdD_fx <- sqrt(diag(VarD_fx))
  
    return(list(Expect=ED_fx,StanDev=sdD_fx))
}
out <- emulator_2d(x_j=x_j,x=x,D=D,beta_0=0,sig2=1,theta=c(1.5,1.5))
ED_fx <- matrix(out$Expect,nrow=length(x1seq))
sdD_fx <- matrix(out$StanDev,nrow=length(x1seq))

library(wesanderson)
#if we don't manually set the number of levels in the contour plot
# as we have done here, then we need to adjust the number of levels
# that we take from the colour palette 
cols <- (wes_palette("Zissou1",19,"continuous"))

#plot the emulator expectation over the 2d space:
filled.contour(x1seq,x2seq,ED_fx,col=cols,xlab=expression(x[1]),ylab=expression(x[2]),
               plot.axes = {axis(1);axis(2);points(x_j,pch=16)},
               lev=seq(min(ED_fx),max(ED_fx),len=20))


#true value of the function over the 2d space: 
real_f <- matrix(2*sin(-0.5*x[,1])+cos(1.5*x[,2])	,nrow=length(x1seq))
#with points
filled.contour(x1seq,x2seq,real_f,col=cols,xlab = expression(x[1]), ylab= expression(x[2]),
               plot.axes = {axis(1);axis(2);points(x_j,pch=16)},lev=seq(min(ED_fx),max(ED_fx),len=20))
#without points
filled.contour(x1seq,x2seq,real_f,col=cols,xlab = expression(x[1]), ylab= expression(x[2]),
               plot.axes = {axis(1);axis(2)},lev=seq(min(ED_fx),max(ED_fx),len=20))

#standard deviation 
cols3 <- (wes_palette("Zissou1",19,"continuous"))
filled.contour(x1seq,x2seq,sdD_fx,col=cols3,xlab = expression(x[1]), ylab= expression(x[2]),
               plot.axes = {axis(1);axis(2);points(x_j,pch=16)})

#Emulator diagnostics example
#To show pinching, run above code with following setup
out <- emulator_2d(x_j=x_j,x=x,D=D,beta_0=0,sig2=1,theta=c(0.5,0.5))

#history matching 1D example
plot_emulator_1d_pastelcolours(ED_fx=out$Expect,sdD_fx=out$StanDev,z=-0.5,sigma_e=0.02,sigma_epsilon=0,yloc=-1.35,ylimit=c(-1.3,1.3),plottype="imp")

#Perform a 2nd wave of runs to improve emulator accuracy at non-implausible points
x_j <- c(1,2,2.2,2.7,3,4)
nx <- 200											
x <- seq(0.7,4.3,len=nx)					
fx_j <- 0.3*exp(((x_j)/2))-1.5				
D <- fx_j	
D
plot(D,pch=16,cex=1.5,ylim=c(-1.3,1.3),lwd=2,xlab="Input Parameter x",ylab="Output f(x)")
out <- emulator_1d(x_j=x_j,x=x,D=D,beta_0=0,sig2=1,theta=1.5)
plot_emulator_1d_pastelcolours(ED_fx=out$Expect,sdD_fx=out$StanDev,z=-0.5,sigma_e=0.02,sigma_epsilon=0,yloc=-1.35,ylimit=c(-1.3,1.3),plottype="emul")
plot_emulator_1d_pastelcolours(ED_fx=out$Expect,sdD_fx=out$StanDev,z=-0.5,sigma_e=0.02,sigma_epsilon=0,yloc=-1.35,ylimit=c(-1.3,1.3),plottype="imp")
plot_emulator_1d_pastelcolours(ED_fx=out$Expect,sdD_fx=out$StanDev,z=-0.5,sigma_e=0.02,sigma_epsilon=0,yloc=-1.4,ylimit=c(-1.3,1.3),plottype="emul_imp")
legend("topleft", legend=c("Adjusted expectation", "3SD credible interval" ),
       col=c("grey40", "grey60"), lty=c(1,1),lwd=c(2,2), cex=0.8)

#history matching example for 2d case
#Define an observation z,observation error sigma_e and model discrepancy sigma_epsilon
z <- -0.5			
sigma_e <- 0.025				# something reasonable e.g. 5% of |z| 
sigma_epsilon <- 0			# zero model discrepancy 

#Calculate the Implausibility 
Imp <- sqrt(  (ED_fx - z)^2 / ( sdD_fx^2 + sigma_epsilon^2 + sigma_e^2)  )
#plot the implausibility
levs <- c(seq(0,5,0.5),max(Imp))						# levels to colour the implausibility plot
zizcols2 <- wes_palette("Zissou1",11,"continuous")
filled.contour(x1seq,x2seq,Imp,col=zizcols2,lev=levs,xlab=expression(x[1]),ylab=expression(x[2]),
               plot.axes = {axis(1);axis(2);points(x_j,pch=16)})

#add in the non-implausible points from x in white (the points having I(x)<3) 
filled.contour(x1seq,x2seq,Imp,col=zizcols2,lev=levs,xlab=expression(x[1]),ylab=expression(x[2]),
               plot.axes = {axis(1);axis(2);points(x_j,pch=16);points(x[Imp<3,],pch=16,col='white')})

#select a set of 15 wave 2 runs "x_j_wave2" by choosing at random from the non-implausible points in x: 
x_j_wave2 <- x[Imp<3,][sample(sum(Imp<3),15),]

#plot the locations of these new proposed runs in yellow
filled.contour(x1seq,x2seq,Imp,col=zizcols2,lev=levs,xlab=expression(x[1]),ylab=expression(x[2]),
               plot.axes = {axis(1);axis(2);points(x_j,pch=16);points(x[Imp<3,],pch=16,col='white');
                 points(x_j_wave2,pch=16,col="black",cex=1.6);points(x_j_wave2,pch=16,col="yellow",cex=1)})

#combine wave 1 and wave 2 inputs together
x_j_wave1_2 <- rbind(x_j,x_j_wave2)

#run the model for the wave 2 runs and combine with wave 1 
fx_j_wave2 <- 2*sin(-0.5*x_j_wave2[,1])+cos(1.5*x_j_wave2[,2])	
D_wave1_2 <- c(fx_j,fx_j_wave2)
out <- emulator_2d(x_j=x_j_wave1_2,x=x,D=D_wave1_2,beta_0=0,sig2=1,theta=c(1.5,1.5))
ED_fx <- matrix(out$Expect,nrow=length(x1seq))
sdD_fx <- matrix(out$StanDev,nrow=length(x1seq))

#plot of new emulator expectation
filled.contour(x1seq,x2seq,ED_fx,col=cols,
               plot.axes = {axis(1);axis(2);points(x_j,pch=16);points(x_j_wave2,pch=16,col="yellow",cex=1)},
               lev=seq(min(ED_fx),max(ED_fx),len=20))

#calculate and plot wave 2 imp
Imp_w2 <- sqrt(  (ED_fx - z)^2 / ( sdD_fx^2 + sigma_epsilon^2 + sigma_e^2)  )
filled.contour(x1seq,x2seq,Imp_w2,col=zizcols2,lev=levs,xlab=expression(x[1]),ylab=expression(x[2]),#main="Wave 2 Implausibility",
               plot.axes = {axis(1);axis(2);points(x_j,pch=16);points(x_j_wave2,pch=16,col="yellow",cex=1)})

#LHS in higher dimensions
xrange <- rbind(c(900,1100),c(500,700),c(400,600))				
nx_j <- 75
x_j <- lhs(n=nx_j,rect=xrange)				#3D points where we will run our computer model f(x_j) from latin hypercube

library(plot3D)
scatter3D(x_j[,1],x_j[,2],x_j[,3],col='#619CFF',pch=19,cex=0.5,bty = "b2",ticktype="detailed")
#can we plot this in ggplot (bit more aesthetically pleasing)
library(ggplot2)
pairsvars<- as.data.frame(x_j)
colnames(pairsvars) <- c("x","y","z")
pairsvars

a<- ggplot(pairsvars,aes(x=x))+
  labs(x='',y='')+
  geom_histogram(aes(y=..density..), colour="black", fill="white")+
  geom_density(alpha=.2, fill="#FF6666") +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),axis.ticks = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
b<-ggplot(pairsvars,aes(x=y))+
  labs(x="",y="")+
  geom_histogram(aes(y=..density..), colour="black", fill="white")+
  geom_density(alpha=.2, fill="#FF6666") +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),axis.ticks = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
c<-ggplot(pairsvars,aes(x=z))+
  labs(y="")+
  geom_histogram(aes(y=..density..), colour="black", fill="white")+
  geom_density(alpha=.2, fill="#FF6666") +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),axis.ticks = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank())

d<-ggplot(pairsvars,aes(x=x,y=y)) + 
  geom_point(colour="#619CFF")+
  
  theme(legend.position = "none",axis.text.x = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(x="")
d
e<-ggplot(pairsvars,aes(x=x,y=z)) + 
  geom_point(colour="#619CFF")+
  theme(legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank())

f<-ggplot(pairsvars,aes(x=y,y=z)) + 
  geom_point(colour="#619CFF")+
  theme(legend.position = "none",
        axis.text.y = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(y="")

#note red colour from https://bio304-class.github.io/bio304-fall2017/ggplot-bivariate.html
library(cowplot)
plot_grid(a,NULL,NULL,d,b,NULL,e,f,c, nrow=3)

#For larger dimensions, using facet grid or something else may be more efficient for plotting similar visuals in GGPLOT. 
#For the remainder of this project, however, cowplot is utilised.

library(ContourFunctions)
library(viridis)
#function to plot
sine3dfct <- function(x) {
  10*sin(pi*x[1]*x[3]) + 20*(x[2]-.5)^2 
}
cf_highdim(sine3dfct, 3, color.palette=viridis)

#5d plot
#very inefficient way to produce this 
#10*cos(2.5*pi*x1*x2)+ 20*(x3-0.5)^2+ 10*x4+ sin(x5)
#note that 0.9 is basically instant, but 0.2 takes about 10 minutes
x1 =  seq(-3.14, 3.14, by = 0.2) #len out 31.4
x1 =  seq(-3.14, 3.14, by = 0.9) #len out 7
#x1=seq(-3.14,3.14,length.out = 20)
fgrid <- expand.grid(x1=x1,x2=x1,x3=x1,x4=x1,x5=x1)
head(fgrid)
fgrid$z1 <- 10*cos(2.5*pi*fgrid$x1*fgrid$x2)
fgrid$z2 <- 10*cos(2.5*pi*fgrid$x1)+ 20*(fgrid$x3-0.5)^2
fgrid$z3 <- 10*cos(2.5*pi*fgrid$x2)+ 20*(fgrid$x3-0.5)^2
fgrid$z4 <- 10*cos(2.5*pi*fgrid$x1)+ 10*fgrid$x4
fgrid$z5 <- 10*cos(2.5*pi*fgrid$x2)+ 10*fgrid$x4
fgrid$z6 <- 20*(fgrid$x3-0.5)^2+ 10*fgrid$x4
fgrid$z7 <- 10*cos(2.5*pi*fgrid$x1)+ sin(fgrid$x5)
fgrid$z8 <- 10*cos(2.5*pi*fgrid$x2)+ sin(fgrid$x5)
fgrid$z9 <- 20*(fgrid$x3-0.5)^2+ sin(fgrid$x5)
fgrid$z10 <- 10*fgrid$x4+ sin(fgrid$x5)

plota<- ggplot(fgrid,aes(x=x1,y=x2,z=z1))+
  geom_contour_filled(aes(fill = stat(level))) +
  #theme_void() +
  labs(x = "")+
  #theme_minimal()+
  theme(legend.position = "none",axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),panel.background = element_rect(fill = "white", colour = "white"))  +ylab(expression(x[2]))
plota

#this plot looks a bit dodge? Compare to filled.contour
a <- expand.grid(seq(-4,4,length.out = 100),seq(-4,4,length.out = 100))
a
b <- matrix(10*cos(2.5*pi*a[,1]*a[,2]), 100)
b
filled.contour(x=seq(-4,4,length.out = 100),y=seq(-4,4,length.out = 100), z = b)
?expand.grid


plotb<- ggplot(fgrid,aes(x=x1,y=x3,z=z2))+
  geom_contour_filled(aes(fill = stat(level))) +
  labs(x = "")+
  #theme_minimal()+
  theme(legend.position = "none",axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),panel.background = element_rect(fill = "white", colour = "white")) +ylab(expression(x[3]))


plotc<- ggplot(fgrid,aes(x=x2,y=x3,z=z3))+
  geom_contour_filled(aes(fill = stat(level))) +
  labs(y = "")+ labs(x = "")+
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),panel.background = element_rect(fill = "white", colour = "white")) 
plotc

plotd<- ggplot(fgrid,aes(x=x1,y=x4,z=z4))+
  geom_contour_filled(aes(fill = stat(level))) +
  labs(x= "")+
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),panel.background = element_rect(fill = "white", colour = "white")) +ylab(expression(x[4]))


plote<- ggplot(fgrid,aes(x=x2,y=x4,z=z5))+
  geom_contour_filled(aes(fill = stat(level))) +
  labs(y = "")+labs(x = "")+
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),panel.background = element_rect(fill = "white", colour = "white"))


plotf<- ggplot(fgrid,aes(x=x3,y=x4,z=z6))+
  geom_contour_filled(aes(fill = stat(level))) +
  labs(y = "")+labs(x = "")+
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),panel.background = element_rect(fill = "white", colour = "white"))



plotg<- ggplot(fgrid,aes(x=x1,y=x5,z=z7))+
  geom_contour_filled(aes(fill = stat(level))) +
  
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),panel.background = element_rect(fill = "white", colour = "white"))+xlab(expression(x[1]))+ylab(expression(x[5]))


ploth<- ggplot(fgrid,aes(x=x2,y=x5,z=z8))+
  geom_contour_filled(aes(fill = stat(level))) +
  labs(y = "")+
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),panel.background = element_rect(fill = "white", colour = "white"))+xlab(expression(x[2]))

ploti<- ggplot(fgrid,aes(x=x3,y=x5,z=z9))+
  geom_contour_filled(aes(fill = stat(level))) +
  labs(y = "")+
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),panel.background = element_rect(fill = "white", colour = "white"))+xlab(expression(x[3]))

plotj<- ggplot(fgrid,aes(x=x4,y=x5,z=z10))+
  geom_contour_filled(aes(fill = stat(level))) +
  labs(y = "")+
  theme(
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),panel.background = element_rect(fill = "white", colour = "white"))+xlab(expression(x[4]))

plotj

plot_grid(plota, NULL,NULL,NULL,
          plotb, plotc,NULL,NULL,
          plotd,plote,plotf,NULL,
          plotg,ploth,ploti,plotj,
          nrow = 4,rel_widths=c(1,1,1,1.3), align = "h")


#3D emulation example
#the friedman plot is over 0,1 so let's try and replicate that
x1seqb <- seq(0,1,len=20)
x2seqb <- seq(0,1,len=20)
x3seqb <- seq(0,1,len=20)
xb <- as.matrix(expand.grid(x1seqb,x2seqb,x3seqb))	
nx_jb <- 200
xrangeb <- rbind(c(0,1),c(0,1),c(0,1))	
x_jb <- lhs(n=nx_jb,rect=xrangeb)	
x_j.dfb <- as.data.frame(x_jb)
fx_jb <- 10*sin(pi*x_jb[,1]*x_jb[,3]) + 20*(x_jb[,2]-.5)^2 	#run the model
Db <- fx_jb

#New 3D emulation function 
emulator_3d <- function(x_j,x,D,beta_0=0,sig2=1,theta=c(1,1,1),delta=10^-5){
  
  ### get numbers of points from supplied x_j and x (which are now arguments of the function)
  nx_j <- nrow(x_j)				# number of points in x_j 
  nx <- nrow(x)					# number of points we want to evaluate the emulator for (apply the BL update for) 
  x_both <- rbind(x,x_j)			# join both x and x_j together (just to make covariance calculation easier)
  n <- nx + nx_j					# total number of old and new points
  
  ### Define prior expectation of f(x) and the vector D = f(x_j) (i.e. before we run the function f(x) at x_j)
  E_fx <- rep(beta_0,nx)			# needed for BL update
  E_D  <- rep(beta_0,nx_j)		# needed for BL update
  
  ### Define prior variances and covariances of f(x) and D = f(x_j)
  Sigma <- sig2 *( (1-delta)* exp( - as.matrix(dist( t( t(x_both)/theta) )^2 ) ) + diag(delta,n))  # bit that makes f(x) a smooth function!
  # So now we should have Sigma a variance matrix, which we cut into three pieces for use in BL update
  Var_fx 	 <- Sigma[1:nx,1:nx]
  Var_D 	 <- Sigma[(nx+1):n,(nx+1):n]
  Cov_fx_D <- Sigma[1:nx,(nx+1):n]
  
  #The Bayes Linear Update equations.
  #They give the 'adjusted' expectation and 'adjusted' variance of f(x) given the run data vector D = f(x_j)
  
  ED_fx <- E_fx + Cov_fx_D %*% solve(Var_D) %*% (D - E_D)						# BL update
  VarD_fx <- Var_fx - Cov_fx_D %*% solve(Var_D) %*% t(Cov_fx_D)				# BL update	
  
  ### Extract the diagonal from full variance matrix and square root as we just want the individual sd's.
  sdD_fx <- sqrt(diag(VarD_fx))
  
  ### return the adjusted expection and just the standard deviations ###
  return(list(Expect=ED_fx,StanDev=sdD_fx))
}



outb <- emulator_3d(x_j=x_jb,x=xb,D=Db,beta_0=0,sig2=10,theta=c(1,1,1)) # use emulator to predict at all the points in the grid 
ED_fxb <- matrix(outb$Expect,nrow=length(x1seqb))
SDD_fxb<- matrix(outb$StanDev,nrow=length(x1seqb))
fctgrid<- expand.grid(x1 =  seq(0, 1, length.out = 20), x2 = seq(0, 1, length.out = 20))#,x3 = seq(400, 600, length.out = 20))
fctgrid$a<- as.vector(avZplot2(20,ED_fxb))
fctgrid$b<- as.vector(avZplot3(20,col_swapper(ED_fxb,20)))
fctgrid$c<- as.vector(avZplot3(20,col_swapper_double(transposer(ED_fxb,20),20)))


#Plot using GGPLOT

fcta<-ggplot()+geom_contour_filled(data=fctgrid,aes(x=x1,y=x2,z=a),breaks=seq(0,15,1))+
  labs(x="",y=expression(x[2]))+
  theme(legend.position = "none",axis.text.x = element_blank(),
        axis.ticks = element_blank(),panel.background = element_rect(fill = "white", colour = "white")) +
  geom_point(data=x_j.dfb, aes(x=V1,y=V2),size=0.3)+theme(axis.title.y = element_text(size=15))

fctb<-ggplot()+geom_contour_filled(data=fctgrid,aes(x=x1,y=x2,z=b),breaks=seq(0,15,1))+
  labs(x=expression(x[1]),y=expression(x[3]))+
  
  theme(legend.position = "none",
        axis.ticks = element_blank(),panel.background = element_rect(fill = "white", colour = "white"))+
  guides(fill = guide_colorsteps(barheight = unit(5, "cm")))+
  geom_point(data=x_j.dfb, aes(x=V1,y=V3),size=0.3)+theme(axis.title.y = element_text(size=15),axis.title.x = element_text(size=15))

fctc<-ggplot()+geom_contour_filled(data=fctgrid,aes(x=x1,y=x2,z=c),breaks=seq(0,15,1))+
  labs(y="",x=expression(x[2]))+
  theme(axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = "white", colour = "white")) +
  guides(fill = guide_colorsteps(barheight = unit(5, "cm"),title=""))+
  geom_point(data=x_j.dfb, aes(x=V2,y=V3),size=0.3)+theme(axis.title.x = element_text(size=15))



plot_grid(fcta,NULL,fctb,fctc,rel_widths = c(1.75, 2))

#############################################
##################################################
#section 4: Modelling liquidity

barcoloursmanual <- c("#81d381","#6b9aea","#ff6961","#6b9aea","#ff6961")
LRcheck <- function(Lat, Inf14, Out14, Inf44, Out44){
  LR <- (Lat+Inf14-Out14)/(Out44-Inf44)
  barplot(c(Lat,Inf14,Out14,Inf44,Out44), 
          names.arg=c('L', expression(I[1:14]), expression(O[1:14]), expression(I[15:44]), expression(O[15:44])), 
          col = barcoloursmanual, border = NA)
  c(LR, Lat+Inf14-Out14, Out44-Inf44)
}
LRcheck(500,325,400,650,1000)


nrep <- function(n, mLat, mInf14, mOut14, mInf44, mOut44, sdperc){ #taking mean values, n reps, sd percentage of mean
  Inf14 = rnorm(n, mInf14, mInf14*sdperc) 
  Out14 = rnorm(n, mOut14, mOut14*sdperc)
  Inf44 = rnorm(n, mInf44, mInf44*sdperc)
  Out44 = rnorm(n, mOut44, mOut44*sdperc)
  LR <- (mLat + Inf14 - Out14) / (Out44 - Inf44)
  LR
}

fig4.2i<- nrep(10,500,325,400,650,1000,0.05)
fig4.2ii<- nrep(100,500,325,400,650,1000,0.05)
fig4.2iii<- nrep(10000,500,325,400,650,1000,0.05)

par(mfrow=c(1,3))
hist(fig4.2i,25, main = NULL, col = 'grey',lty="blank", xlab = expression(gamma))
hist(fig4.2ii,55,main = NULL, col = 'grey',lty="blank",xlab = expression(gamma))
hist(fig4.2iii,55,main = NULL, col = 'grey',lty="blank",xlab = expression(gamma))



Out44 = as.data.frame(rnorm(10000, 1000, 1000*0.05))
colnames(Out44)= "Out44"
pa<-ggplot(Out44,aes(x=Out44))+
  geom_boxplot(outlier.size = 0.04,outlier.colour = 'red',outlier.shape = 8)+
  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(x=expression(O[15:44]))

pb<-ggplot(Out44,aes(sample=Out44))+
  stat_qq()+
  theme(
    
    panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(x="Theoretical quantiles",y="Sample quantiles")

aaa<-(as.data.frame(fig4.2iii))
pc<-ggplot(aaa,aes(x=fig4.2iii))+
  geom_boxplot(outlier.size = 0.04,outlier.colour = 'red',outlier.shape = 8)+
  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(x=expression(gamma))


pd<-ggplot(aaa,aes(sample=fig4.2iii))+
  stat_qq()+
  theme(
    
    panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(x="Theoretical quantiles",y="Sample quantiles")
plot_grid(pa,NULL,pb,pc,NULL,pd,nrow = 2, rel_widths = c(1, 0.1, 1))


#Cauchy
par(mfrow=c(1,1))
curve(dcauchy(x,0, 0.5),from = -6, to = 6, col = "coral", xlab = "x", ylab = "Density", lwd = 2)
curve(dcauchy(x,0, 1),from = -6, to = 6, , col = "grey", add = TRUE, lwd = 2)
#curve(dcauchy(x,0, 2),from = -6, to = 6, col = " green", add = TRUE, lwd = 2)
curve(dcauchy(x,-2, 1),from = -6, to = 6, col = "blue", add = TRUE, lwd = 2)
curve(dcauchy(x,-1, 0.75),from = -6, to = 6, col = "green", add = TRUE, lwd = 2)
legend("topleft", legend = paste0("α = ", c(0, 0, -1, -2), ", β=", c(0.5,1,0.75,1)),
       col = c("coral", "grey ", " green", "blue"),
       lty = 1, lwd = 2)


nrepmean <- function(N, n, mLat, mInf14, mOut14, mInf44, mOut44, sdperc){
  out <- matrix( ncol=N, nrow=1)
  for (i in 1:N){
    out[,i] <- mean(nrep(n, mLat, mInf14, mOut14, mInf44, mOut44, sdperc))
  }
  out
}

par(mfrow=c(1,3))
figure4.5i<- nrepmean(10, 25, 500, 325, 400, 650, 1000, 0.05)
figure4.5ii<- nrepmean(100, 25, 500, 325, 400, 650, 1000, 0.05)
figure4.5iii<- nrepmean(10000, 25, 500, 325, 400, 650, 1000, 0.05)
hist(figure4.5i,25, main = NULL, col = 'grey',lty="blank", xlab = expression(gamma))
hist(figure4.5ii,55, main = NULL, col = 'grey',lty="blank", xlab = expression(gamma))
hist(figure4.5iii,25, main = NULL, col = 'grey',lty="blank", xlab = expression(gamma))

library(MASS)
library(psych)
Mu <- c(0,0,0)
Sig <- matrix(c(1, -0.7, -.5,
                -0.7, 1, 0.6,
                -0.5, 0.6, 1), 
              nrow=3)
X <- mvrnorm(1000, mu=Mu, Sigma = Sig, 
             empirical = TRUE)  
pairs.panels(X, hist.col = "grey")

#Plot this ourselves using GGPLOT
mvn1 <- as.data.frame(mvrnorm(150, mu = Mu, Sigma = Sig))
colnames(mvn1) <- c("x1","x2","x3")
mvn1

a<- ggplot(mvn1,aes(x=x1))+
  labs(x='',y='')+
  geom_histogram(aes(y=..density..), colour="black", fill="white")+
  geom_density(alpha=.2, fill="#FF6666") +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),axis.ticks = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
b<-ggplot(mvn1,aes(x=x2))+
  labs(x="",y="")+
  geom_histogram(aes(y=..density..), colour="black", fill="white")+
  geom_density(alpha=.2, fill="#FF6666") +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),axis.ticks = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
c<-ggplot(mvn1,aes(x=x3))+
  labs(x=expression(x[3]),y="")+
  geom_histogram(aes(y=..density..), colour="black", fill="white")+
  geom_density(alpha=.2, fill="#FF6666") +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),axis.ticks = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank())

d<-ggplot(mvn1,aes(x=x1,y=x2)) + 
  stat_density_2d(aes(fill = ..level..), geom = "polygon") +
  scale_fill_continuous(low="lavenderblush", high="red") +
  geom_jitter(alpha=0.5, size = 1.1) +
  theme(legend.position = "none",axis.text.x = element_blank(),
        axis.text.y = element_blank(),axis.ticks = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(x="",y=expression(x[2]))

e<-ggplot(mvn1,aes(x=x1,y=x3)) + 
  stat_density_2d(aes(fill = ..level..), geom = "polygon") +
  scale_fill_continuous(low="lavenderblush", high="red") +
  geom_jitter(alpha=0.5, size = 1.1) +
  theme(legend.position = "none",axis.text.x = element_blank(),
        axis.text.y = element_blank(),axis.ticks = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(x=expression(x[1]),y=expression(x[3]))

f<-ggplot(mvn1,aes(x=x2,y=x3)) + 
  stat_density_2d(aes(fill = ..level..), geom = "polygon") +
  scale_fill_continuous(low="lavenderblush", high="red") +
  geom_jitter(alpha=0.5, size = 1.1) +
  theme(legend.position = "none",axis.text.x = element_blank(),
        axis.text.y = element_blank(),axis.ticks = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(x=expression(x[2]),y="")
#note colour scheme from https://bio304-class.github.io/bio304-fall2017/ggplot-bivariate.html


plot_grid(a,NULL,NULL,d,b,NULL,e,f,c, nrow=3)


#######



#gaussian copula

bvn <- function(n,mu1,mu2,sd1,sd2,r){
  mu <- c(mu1,mu2)
  sigma <- matrix(c(sd1^2, sd1*sd2*r, sd1*sd2*r, sd2^2), nrow=2)
  A <- t(chol(sigma))
  Z <- matrix(rnorm(2*n),2,n) #just N(1,0)
  bvnout <- t(A %*% Z) + matrix(rep(mu,n), byrow=TRUE,ncol=2)
  bvnout
}

copula<- pnorm(bvn(1000,0,0,1,1,0.7))
par(mfrow=c(1,3))
plot(copula)
copula[,1]
hist(copula[,1],col='royal blue',freq = FALSE,50)
abline(h=1,col='red')
hist(copula[,2],col='royal blue',freq = FALSE,50)
abline(h=1,col='red')



#hist(-2+tan(pi*(copula[,1]-0.5)),20)
par(mfrow=c(1,2))
hist(qnorm(copula[,2],600,30),100,col='royal blue',freq = FALSE)
curve(dnorm(x, 600, 30),from = 300, to = 700, col = "red", add = TRUE, lwd = 2)


cauchydata<-(qcauchy(copula[,1],800,20))
#hist(cauchydata,200)
hist(cauchydata,100,col='royal blue',freq = FALSE)
curve(dcauchy(x, 800, 20),from = 0, to = 2000, col = "red", add = TRUE, lwd = 2)

#plot these in GGPLOT
copula.df <- as.data.frame(copula)
fig4.8i<-ggplot(copula.df,aes(x=V1,y=V2)) + 
  geom_point(colour="#619CFF")+
  
  theme(legend.position = "none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(x=expression(U[1]),y=expression(U[2]))+theme_classic()

fig4.8ii <-   ggplot(copula.df,aes(x=V1)) + 
  geom_histogram(aes(y=..density..), colour="black", fill="#619CFF")+theme_classic()+
  geom_hline(yintercept=1, linetype="dashed", color = "#F8766D")+ labs(x=expression(U[1]),y="Density")
  
fig4.8iii <-   ggplot(copula.df,aes(x=V2)) + 
  geom_histogram(aes(y=..density..), colour="black", fill="#619CFF")+theme_classic()+
  geom_hline(yintercept=1, linetype="dashed", color = "#F8766D")+ labs(x=expression(U[2]),y="Density")

plot_grid(fig4.8i,fig4.8ii,fig4.8iii,nrow=1)




copula.df$cauchydata <- (qcauchy(copula[,1],800,20))
cauchydata.df<- as.data.frame(cauchydata)
fig4.9i<-ggplot(cauchydata.df,aes(x=cauchydata)) + 
  geom_histogram(aes(y=..density..),binwidth = 30, colour="black", fill="#619CFF")+theme_classic()+
  stat_function(fun = function(x) dcauchy(x, 800, 20) ,
                color = "#F8766D", size = 1)+ labs(x=expression(X[1]),y="Density")
  #geom_hline(yintercept=1, linetype="dashed", color = "#F8766D")+ labs(x=expression(U[1]),y="Density")
copula.df$X2<- qnorm(copula.df$V2,600,30)
fig4.9ii<-ggplot(copula.df,aes(x=X2)) + 
  geom_histogram(aes(y=..density..),binwidth = 8, colour="black", fill="#619CFF")+theme_classic()+
  stat_function(fun = function(x) dnorm(x, 600, 30) ,
                color = "#F8766D", size = 1)+labs(x=expression(X[2]),y="Density")



plot_grid(fig4.9i,fig4.9ii,nrow=1)



newLRcheck <- function(x_1,x_2,x_3,x_4,x_5,x_6,x_7){
  LR<-(x_1+x_2-x_3)/(x_4+x_5+x_6-x_7)
  barplot(c(x_1,x_2,x_3,x_4,x_5,x_6,x_7), names.arg=c(expression(x[1]),expression(x[2]),expression(x[3]),expression(x[4]),expression(x[5]),expression(x[6]),expression(x[7])))
  a= c(LR, x_1+x_2-x_3, x_4+x_5+x_6-x_7)
  a
}

newLRcheck(1000,3000,2000,2000,1200,800,2340)

ratio_corr <- function(n,x_1,mx_2,mx_3,mx_4,mx_5,x_6,x_7,sdperc,rx_3x_4){
  x_2=rnorm(n,mx_2,mx_2*sdperc)
  x_3=bvn(n,mx_3,mx_4,mx_3*sdperc,mx_4*sdperc,rx_3x_4)[,1]
  x_4=bvn(n,mx_3,mx_4,mx_3*sdperc,mx_4*sdperc,rx_3x_4)[,2]
  x_5=rnorm(n,mx_5,mx_5*sdperc)
  LR<-(x_1+x_2-x_3)/(x_4+x_5+x_6-x_7)
  LR
}

m<-ratio_corr(1000,1000,3000,2000,2000,1200,800,2340,00.01,0.01)
hist(m)
m.df<- (as.data.frame(m))
histogram7dimplot<- ggplot(m.df,aes(x=m))  +geom_histogram(aes(y=..density..),colour="black", fill="#619CFF")+
  theme_classic()+labs(x=expression(gamma),y="Density")
  


###############################################
#Section 5
#managing inflows to reach a target ratio
invrep <- function(n, LR, mLat, mOut14, mInf44, mOut44, sdperc){ #taking mean values, n reps, sd percentage of mean
  Out14 = rnorm(n, mOut14, mOut14*sdperc)
  Inf44 = rnorm(n, mInf44, mInf44*sdperc)
  Out44 = rnorm(n, mOut44, mOut44*sdperc)
  Inf14 <- LR*(Out44-Inf44) - mLat + Out14
  Inf14
}
invrepmean <- function(N, n, LR, mLat, mOut14, mInf44, mOut44, sdperc){
  out <- matrix( ncol=N, nrow=2)
  for (i in 1:N){
    out[1,i] <- mean(invrep(n, LR, mLat, mOut14, mInf44, mOut44, sdperc))
    out[2,i] <- sd(invrep(n, LR, mLat, mOut14, mInf44, mOut44, sdperc))
  }
  out
  
}
invrepcrit <- function(LR=1, mLat, mOut14, mInf44, mOut44){
  Inf14 <- LR*(mOut44-mInf44) - mLat + mOut14
  Inf14
}

#Note that if we wish to produce the plot with black line and dotted lines for Inf14 actual then uncomment the 
# commented geom_hline functions in lrdist

lrdist <- function(n, mLR, mLat, mOut14, mInf44, mOut44, sdperc,k){ #taking mean values, n reps, sd percentage of mean
  # K the variability the bank can change Inf14 by
  #  as a percentage of the Inf value
  LR = runif(n, mLR-0.05*mLR, mLR+0.05*mLR )
  Out14 = rnorm(n, mOut14, mOut14*sdperc)
  Inf44 = rnorm(n, mInf44, mInf44*sdperc)
  Out44 = rnorm(n, mOut44, mOut44*sdperc)
  Inf14req <- LR*(Out44-Inf44) - mLat + Out14 # this value will get our LR exactly
  Inf14crit <- matrix(ncol=n, nrow=1) # this is the minimum required Inf14 for LR to be >= 100%
  for (i in 1:n){
    Inf14crit[,i] <- invrepcrit(1,mLat,Out14[i],Inf44[i],Out44[i])
  }
  x <- c(1:n)
  InfSort = sort(Inf14req,index.return=TRUE)
  Inf14act <- mLR*(mOut44-mInf44) - mLat + mOut14 
  print(Inf14act)
    a<-as.data.frame(cbind(x,Inf14req[InfSort$ix],Inf14crit[InfSort$ix]))
  colnames(a)<- c("x","infreq","infcrit")
  b<- ggplot()+geom_line(aes(x=x,y=infreq),colour="#81d381",a)+ geom_line(aes(x=x, y=infcrit, colour='blue'), a)+
    #geom_hline(yintercept=Inf14act, color = "black")+
    #geom_hline(yintercept=Inf14act*(1-k), linetype="dashed", color = "black")+
    #geom_hline(yintercept=Inf14act*(1+k), linetype="dashed", color = "black")+
    ylab(TeX("$I_{1:14}$"))+xlab("Number of simulations")+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position = "")
  return(b)
  
}

lrdist(10000,1.2, 500, 400, 650, 1000, 0.05,0.1)
#add a bit more sd
lrdist(10000,1.2, 500, 400, 650, 1000, 0.15,0.1)
# a bit less sd
lrdist(10000,1.2, 500, 400, 650, 1000, 0.025,0.1)

#Function to do trace plot!
lrdistb <- function(n, mLR, mLat, mOut14, mInf44, mOut44, sdperc,k){ #taking mean values, n reps, sd percentage of mean
  # K the variability the bank can change Inf14 by
  #  as a percentage of the Inf value
  LR = runif(n, mLR-0.05*mLR, mLR+0.05*mLR )
  Out14 = rnorm(n, mOut14, mOut14*sdperc)
  Inf44 = rnorm(n, mInf44, mInf44*sdperc)
  Out44 = rnorm(n, mOut44, mOut44*sdperc)
  Inf14req <- LR*(Out44-Inf44) - mLat + Out14 # this value will get our LR exactly
  Inf14crit <- matrix(ncol=n, nrow=1) # this is the minimum required Inf14 for LR to be >= 100%
  for (i in 1:n){
    Inf14crit[,i] <- invrepcrit(1,mLat,Out14[i],Inf44[i],Out44[i])
  }
  x <- c(1:n)
  InfSort = sort(Inf14req,index.return=TRUE)
  Inf14act <- mLR*(mOut44-mInf44) - mLat + mOut14 
  print(Inf14act)
  a<-as.data.frame(cbind(x,Inf14req[x],Inf14crit[x]))
  colnames(a)<- c("x","infreq","infcrit")
  b<- ggplot()+geom_line(aes(x=x,y=infreq),colour="#81d381",a)+ geom_line(aes(x=x, y=infcrit, colour='blue'), a)+
    ylab(TeX("$I_{1:14}$"))+xlab("Number of simulations")+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position = "")
  return(b)
  
}
lrdistb(10000,1.2, 500, 400, 650, 1000, 0.05,0.1)

#Function to show proportions of failures to hit target
lrdist2 <- function(n, mLR, mLat, mOut14, mInf44, mOut44, sdperc,k){ #taking mean values, n reps, sd percentage of mean
  # K the variability the bank can change Inf14 by
  LR = runif(n, mLR-0.05*mLR, mLR+0.05*mLR )
  Out14 = rnorm(n, mOut14, mOut14*sdperc)
  Inf44 = rnorm(n, mInf44, mInf44*sdperc)
  Out44 = rnorm(n, mOut44, mOut44*sdperc)
  Inf14req <- LR*(Out44-Inf44) - mLat + Out14
  Inf14crit <- matrix(ncol=n, nrow=1)
  InfSort = sort(Inf14req,index.return=TRUE)
  
  for (i in 1:n){
    Inf14crit[,i] <- invrepcrit(1,mLat,Out14[i],Inf44[i],Out44[i])
  }
  x <- c(1:n)
  Inf14act <- mLR*(mOut44-mInf44) - mLat + mOut14 
  greenv <- matrix(ncol=n,nrow=1)
  bluev <- matrix(ncol=n,nrow=1)
  greenv2 = as.numeric(Inf14req>Inf14act*(1+k))
  bluev = as.numeric(Inf14crit>Inf14act*(1+k))
  totg <- sum(greenv) #this is how many times we cannot hit LR = 120%
  totb <- sum(bluev) #this is how many times we cannot hit LR = 100%
  print(c(totb,totb/n))
  #totb
}
lrdist2(250,1.2, 500, 400, 650, 1000, 0.05,0.1)
lrdist2(1000,1.2, 500, 400, 650, 1000, 0.05,0.1)
#lrdist2(10000,1.2, 500, 400, 650, 1000, 0.05,0.3)
lrdist2(100000,1.2, 500, 400, 650, 1000, 0.05,0.1)
lrdist2(2000000,1.2, 500, 400, 650, 1000, 0.05,0.1)
#lrdist2(100000,1.2, 1000, 2000, 2340, 4000,0.05,0.1)
#lrdist2(200000,1.2, 1000, 2000, 2340, 4000,0.05,0.1)
a2

#Function to explore quantiles of the function

lrdist3 <- function(n, mLR, mLat, mOut14, mInf44, mOut44, sdperc,k,percentage){ 
  #taking mean values, n reps, sd percentage of mean
  #note that the percentage should be inputted as a vector
  # K the variability the bank can change Inf14 by
  LR = runif(n, mLR-0.05*mLR, mLR+0.05*mLR )
  Out14 = rnorm(n, mOut14, mOut14*sdperc)
  Inf44 = rnorm(n, mInf44, mInf44*sdperc)
  Out44 = rnorm(n, mOut44, mOut44*sdperc)
  Inf14req <- LR*(Out44-Inf44) - mLat + Out14
  Inf14crit <- matrix(ncol=n, nrow=1)
  InfSort = sort(Inf14req,index.return=TRUE)
  
  for (i in 1:n){
    Inf14crit[,i] <- invrepcrit(1,mLat,Out14[i],Inf44[i],Out44[i])
  }
  x <- c(1:n)
  bluev <- matrix(ncol=n,nrow=1)
  bluev<- quantile(as.numeric(Inf14crit),probs = (percentage))
  print(bluev)
  
}
lrdist3(10000,1.2, 500, 400, 650, 1000, 0.05,0.1,c(0.9995,0.999,0.995,0.99,0.9))
#  99.95%    99.9%    99.5%      99%      90% 
#455.6443 440.4023 416.6786 400.3790 332.3886 


#####

#Emulating liquidity models
#1D Emulation example

#Now let's try to emulate our liquidity forecasting model
#We use 1D emulation code as in Chapter 2: refer there for example code
emulator_1d <- function(x_j,x,D,beta_0=0,sig2=1,theta=15){
  
  nx_j <- length(x_j)			
  nx <- length(x)					
  x_both <- c(x,x_j)			
  n <- nx + nx_j					
  
  E_fx <- rep(beta_0,nx)			
  E_D <- rep(beta_0,nx_j)			
  
  Sigma <- sig2 * exp( - as.matrix(dist(x_both))^2/theta^2 )   
  
  Var_fx 	 <- Sigma[1:nx,1:nx]
  Var_D 	 <- Sigma[(nx+1):n,(nx+1):n]
  Cov_fx_D <- Sigma[1:nx,(nx+1):n]							
  
  ED_fx <- E_fx + Cov_fx_D %*% solve(Var_D) %*% (D - E_D)							
  VarD_fx <- Var_fx - Cov_fx_D %*% solve(Var_D) %*% t(Cov_fx_D)					
  sdD_fx <- sqrt(diag(VarD_fx))
  return(list(Expect=ED_fx,StanDev=sdD_fx))
}
plot_emulator_1d_pastelcolours <- function(ED_fx,sdD_fx,z=-0.5,sigma_e=0.02,sigma_epsilon=0,yloc=-1.35,ylimit,plottype="emul_imp"){
  
  if(plottype=="emul"){
    plot(x,ED_fx,ty="l",col="#6b9aea",ylim=ylimit,lwd=2,xlab=expression(I[i:j]),ylab="Emulator Output")
    lines(x,ED_fx+2*sdD_fx,col="#ff6961",lwd=2)
    lines(x,ED_fx-2*sdD_fx,col="#ff6961",lwd=2)
    points(x_j,D,pch=16,cex=1.5)				
  }
  
  if(plottype=="imp"){
    Imp <- sqrt(  (ED_fx - z)^2 / ( sdD_fx^2 + sigma_epsilon^2 + sigma_e^2)  )
    plot(x,Imp,ty="l",ylim=c(0,30),xlab="Input Parameter x",ylab="Implausibility",lwd=2)
    abline(h=3,col="#6b9aea",lwd=2)
    
    points(x[Imp>5],rep(yloc,sum(Imp>5)),pch=16,col="#ff6961",cex=1.5)				
    points(x[(3<Imp) & (Imp<5)],rep(yloc,sum((3<Imp) & (Imp<5))),pch=16,col="#ffb347",cex=1.5)
    points(x[Imp<3],rep(yloc,sum(Imp<3)),pch=16,col="#81d381",cex=1.5)
  }
  
  if(plottype=="emul_imp"){
    
    plot(x,ED_fx,ty="l",col="grey40",ylim=ylimit,lwd=2,xlab="Input Parameter x",ylab="Emulator Output")
    polygon(c(min(x)-1,min(x)-1,max(x)+1,max(x)+1),c(z-3*sigma_e,z+3*sigma_e,z+3*sigma_e,z-3*sigma_e),col="#ffe6e5",border = NA)
    lines(x,ED_fx,col="grey40",lwd=2)
    lines(x,ED_fx+2*sdD_fx,col="grey60",lwd=2)
    
    lines(x,ED_fx-2*sdD_fx,col="gray60",lwd=2)
    points(x_j,D,pch=16,cex=1.5)				
    
    Imp <- sqrt(  (ED_fx - z)^2 / ( sdD_fx^2 + sigma_epsilon^2 + sigma_e^2)  )
    
    
    abline(h=c(z,z+3*sigma_e,z-3*sigma_e),lty=c(1,2,2),lwd=1.4)
    
    
    points(x[Imp>5],rep(yloc,sum(Imp>5)),pch=16,col="#ff6961",cex=1.5)				
    points(x[(3<Imp) & (Imp<5)],rep(yloc,sum((3<Imp) & (Imp<5))),pch=16,col="#ffb347",cex=1.5)
    points(x[Imp<3],rep(yloc,sum(Imp<3)),pch=16,col="#81d381",cex=1.5)
    print((x[Imp<3]))
  }
}

####                                          LAT
#Consider simple lcr equation  gamma =  ---------------
#                                       Outflows+Inflows    
#keep LAT const, outflows = const, inflows = x



#2D example
#Follow the method of Chapter 2 using the emulator_2d function and wesanderson colour palette
emulator_2d <- function(x_j,x,D,beta_0=0,sig2=1,theta=c(1,1),delta=10^-5){
  
  ### get numbers of points from supplied x_j and x (which are now arguments of the function)
  nx_j <- nrow(x_j)				# number of points in x_j 
  nx <- nrow(x)					# number of points we want to evaluate the emulator for (apply the BL update for) 
  x_both <- rbind(x,x_j)			# join both x and x_j together (just to make covariance calculation easier)
  n <- nx + nx_j					# total number of old and new points
  
  ### Define prior expectation of f(x) and the vector D = f(x_j) (i.e. before we run the function f(x) at x_j)
  E_fx <- rep(beta_0,nx)			# needed for BL update
  E_D  <- rep(beta_0,nx_j)		# needed for BL update
  
  ### Define prior variances and covariances of f(x) and D = f(x_j)
  Sigma <- sig2 *( (1-delta)* exp( - as.matrix(dist( t( t(x_both)/theta) )^2 ) ) + diag(delta,n))  # bit that makes f(x) a smooth function!
  # So now we should have Sigma a variance matrix, which we cut into three pieces for use in BL update
  Var_fx 	 <- Sigma[1:nx,1:nx]
  Var_D 	 <- Sigma[(nx+1):n,(nx+1):n]
  Cov_fx_D <- Sigma[1:nx,(nx+1):n]
  
  ### The Bayes Linear Update equations.
  ### They give the 'adjusted' expectation and 'adjusted' variance of f(x) given the run data vector D = f(x_j)
  ### (Note: to multiply matrices in R use %*% and to invert a matrix use solve().)
  ED_fx <- E_fx + Cov_fx_D %*% solve(Var_D) %*% (D - E_D)						# BL update
  VarD_fx <- Var_fx - Cov_fx_D %*% solve(Var_D) %*% t(Cov_fx_D)				# BL update	
  
  ### Extract the diagonal from full variance matrix and square root as we just want the individual sd's.
  sdD_fx <- sqrt(diag(VarD_fx))
  
  ### return the adjusted expection and just the standard deviations ###
  return(list(Expect=ED_fx,StanDev=sdD_fx))
}
library(tgp)		#for lhs function
library(wesanderson)

#Note that for 4+ waves of the history matching process, we reapply the same method for carrying out the first wave



#3d emulation example
emulator_3d <- function(x_j,x,D,beta_0=0,sig2=1,theta=c(1,1,1),delta=10^-5){
  
  ### get numbers of points from supplied x_j and x (which are now arguments of the function)
  nx_j <- nrow(x_j)				# number of points in x_j 
  nx <- nrow(x)					# number of points we want to evaluate the emulator for (apply the BL update for) 
  x_both <- rbind(x,x_j)			# join both x and x_j together (just to make covariance calculation easier)
  n <- nx + nx_j					# total number of old and new points
  
  ### Define prior expectation of f(x) and the vector D = f(x_j) (i.e. before we run the function f(x) at x_j)
  E_fx <- rep(beta_0,nx)			# needed for BL update
  E_D  <- rep(beta_0,nx_j)		# needed for BL update
  
  ### Define prior variances and covariances of f(x) and D = f(x_j)
  Sigma <- sig2 *( (1-delta)* exp( - as.matrix(dist( t( t(x_both)/theta) )^2 ) ) + diag(delta,n))  # bit that makes f(x) a smooth function!
  # So now we should have Sigma a variance matrix, which we cut into three pieces for use in BL update
  Var_fx 	 <- Sigma[1:nx,1:nx]
  Var_D 	 <- Sigma[(nx+1):n,(nx+1):n]
  Cov_fx_D <- Sigma[1:nx,(nx+1):n]
  
  #The Bayes Linear Update equations.
  #They give the 'adjusted' expectation and 'adjusted' variance of f(x) given the run data vector D = f(x_j)
  
  ED_fx <- E_fx + Cov_fx_D %*% solve(Var_D) %*% (D - E_D)						# BL update
  VarD_fx <- Var_fx - Cov_fx_D %*% solve(Var_D) %*% t(Cov_fx_D)				# BL update	
  
  ### Extract the diagonal from full variance matrix and square root as we just want the individual sd's.
  sdD_fx <- sqrt(diag(VarD_fx))
  
  ### return the adjusted expection and just the standard deviations ###
  return(list(Expect=ED_fx,StanDev=sdD_fx))
}

#carry out 3d emulation similar to before, but to manipulate our output to plot using GGPLOT
# we use a combination of the following functions


mat_split <- function(M, r, c){
  nr <- ceiling(nrow(M)/r)
  nc <- ceiling(ncol(M)/c)
  newM <- matrix(NA, nr*r, nc*c)
  newM[1:nrow(M), 1:ncol(M)] <- M
  
  div_k <- kronecker(matrix(seq_len(nr*nc), nr, byrow = TRUE), matrix(1, r, c))
  matlist <- split(newM, div_k)
  N <- length(matlist)
  mats <- unlist(matlist)
  dim(mats)<-c(r, c, N)
  return(mats)
} # allows us to split our EDfx matrix over levels of nth variable. 
#We can then compute the average over the nth variable 
#https://stackoverflow.com/questions/11641701/sum-a-list-of-matrices
avZplot2 <- function(k, input_matrix){
  #k is the size of the k^n grid
  a<-mat_split(input_matrix,k,k)
  b<-matrix(0,k,k)
  for(k in 1:k){b<- b+ a[,,k]}
  return(b/k)
}
transposer <- function(matrix,k){
  a<- mat_split(matrix,k,k)
  b<- mat_split(matrix(0,k,k^2),k,k)
  for(k in 1:k){b[,,k]<- t(a[,,k])}
  
  return(b)
}
col_swapper<- function(matrix,k){
  a<- mat_split(matrix,k,k)
  b<- mat_split(matrix(0,k,k^2),k,k)
  for(i in 1:k)for(j in 1:k){b[,i,j]<- a[,j,i]}
  return(b)
}
col_swapper_double<- function(double,k){
  a<- double
  b<- mat_split(matrix(0,k,k^2),k,k)
  for(i in 1:k)for(j in 1:k){b[,i,j]<- a[,j,i]}
  return(b)
  
}
avZplot3 <- function(k, input_double){
  #k is the size of the k^n grid
  a<-input_double
  b<-matrix(0,k,k)
  for(k in 1:k){b<- b+ a[,,k]}
  return(b/k)
}

#We have now made the functions we need to carry out a plot of a 3d emulation
#so we keep L as const 500, Out14 const 400
#predictors are Inf14, Inf44 and Out44
#inf14 N(325,16.25)
#Inf44 N(650,32.5)
#out44 N(1000,50)

#going to attempt to avoid the LHS function and use actual normal distributions for the x_j
#set up the run locations. 3d points where we will run our computer model f(x_j) 
x_jb <- cbind(rnorm(20,325,16.25),rnorm(20,650,32.5),rnorm(20,1000,50))	
x_jb.df <- as.data.frame(x_jb)
#make a grid of points in 3d to evaluate the emulator at
x1seqb <- seq(275,375,len=20)
x2seqb<- seq(550,750,len=20)
x3seqb <- seq(800,1200,len=20)
xb <- as.matrix(expand.grid(x1seqb,x2seqb,x3seqb))	
											
#run the model 
fx_jb <- (500+x_jb[,1]-400)/(x_jb[,3]-x_jb[,2])
mean(fx_jb)
Db <- fx_jb
outb <- emulator_3d(x_j=x_jb,x=xb,D=Db,beta_0=0,sig2=30,theta=c(250,250,250)) # use emulator to predict at all the points in the grid 
ED_fxb <- matrix(outb$Expect,nrow=length(x1seqb))
sdD_fxb <- matrix(outb$StanDev,nrow=length(x1seqb))
mean(ED_fxb)

#Plot the 3 bits of data
EDf_InfInf<-avZplot2(20,ED_fxb)
EDf_inf14out44<-avZplot3(20,col_swapper(ED_fxb,20))
EDf_inf44out44<-avZplot3(20,col_swapper_double(transposer(ED_fxb,20),20))

jose1<- expand.grid(x1 =  seq(275, 375, length.out = 20), x2 = seq(550, 750, length.out = 20))#,x3 = seq(400, 600, length.out = 20))
jose2<- expand.grid(x1 =  seq(275, 375, length.out = 20),x3 = seq(800, 1200, length.out = 20))
jose3<- expand.grid(x2 = seq(550, 750, length.out = 20),x3 = seq(800, 1200, length.out = 20))

jose1$a <- as.vector(EDf_InfInf)
jose2$a <- as.vector(EDf_inf14out44)
jose3$a <- as.vector(EDf_inf44out44)
library(latex2exp)
a<-ggplot()+geom_contour_filled(data=jose1,aes(x=x1,y=x2,z=a),breaks=seq(0.3,2.5,0.1))+geom_point(data=x_jb.df, aes(x=V1,y=V2),size=0.3)+ 
  labs(x="",y=TeX("$I_{15:44}$"))+
  theme(legend.position = "",axis.text.x=element_blank(),panel.background = element_rect(fill = "white", colour = "white"))
b<-ggplot()+geom_contour_filled(data=jose2,aes(x=x1,y=x3,z=a),breaks=seq(0.3,2.5,0.1))+  geom_point(data=x_jb.df, aes(x=V1,y=V3),size=0.3)+ 
  labs(x=TeX("$I_{1:14}$"), y=TeX("$O_{15:44}$"))+
  theme(legend.position = "",panel.background = element_rect(fill = "white", colour = "white"))
c<-ggplot()+geom_contour_filled(data=jose3,aes(x=x2,y=x3,z=a),breaks=seq(0.3,3.5,0.1))+
  geom_point(data=x_jb.df, aes(x=V2,y=V3),size=0.3)+ labs(x=TeX("$I_{15:44}$"), y="")+
  theme(axis.text.y=element_blank(),panel.background = element_rect(fill = "white", colour = "white"))
library(cowplot)
plot_grid(a,NULL,b,c,nrow = 2, rel_widths = c(1,1.48))

#Now proceed as before to carry out implausibility plots and further waves of history matching!



