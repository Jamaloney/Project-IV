#2D example
library(tgp) #for lhs function
library(wesanderson) # for colour palette

#############
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
emulator_2d <- function(x_j,x,D,beta_0=0,sig2=1,theta=c(1,1),delta=10^-5){
  
  
  nx_j <- nrow(x_j)				# number of points in x_j 
  nx <- nrow(x)					# number of points we want to evaluate the emulator for (apply the BL update for) 
  x_both <- rbind(x,x_j)			# join both x and x_j together (just to make covariance calculation easier)
  n <- nx + nx_j					# total number of old and new points
  
  #Define prior expectation of f(x) and the vector D = f(x_j) 
  E_fx <- rep(beta_0,nx)			
  E_D  <- rep(beta_0,nx_j)		
  
  #Prior variances and covariances of f(x) and D = f(x_j)
  Sigma <- sig2 *( (1-delta)* exp( - as.matrix(dist( t( t(x_both)/theta) )^2 ) ) + diag(delta,n))  
  Var_fx 	 <- Sigma[1:nx,1:nx]
  Var_D 	 <- Sigma[(nx+1):n,(nx+1):n]
  Cov_fx_D <- Sigma[1:nx,(nx+1):n]
  
  #Bayes linear update equations.
  ED_fx <- E_fx + Cov_fx_D %*% solve(Var_D) %*% (D - E_D)						
  VarD_fx <- Var_fx - Cov_fx_D %*% solve(Var_D) %*% t(Cov_fx_D)			
  sdD_fx <- sqrt(diag(VarD_fx))
  return(list(Expect=ED_fx,StanDev=sdD_fx))
}

############

#Regression model vs 1D emulation example
#Introduce function to plot our emulation
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


####
#2d emulation of LCR 


xrange <- rbind(c(900,1100),c(500,700))				#range of input space to explore
nx_j <- 75
x_j <- lhs(n=nx_j,rect=xrange)				#2d points where we will run our computer model f(x_j) from latin hypercube
plot(x_j,pch=16)						


#grid of points in 2d to evaluate the emulator at 
x1seq <- seq(900,1100,len=50)
x2seq <- seq(500,700,len=50)
x <- as.matrix(expand.grid(x1seq,x2seq))	

#run the model 
fx_j <- 500/(x_j[,1]-x_j[,2])	
D <- fx_j								
out <- emulator_2d(x_j=x_j,x=x,D=D,beta_0=0,sig2=30,theta=c(55,55))
ED_fx <- matrix(out$Expect,nrow=length(x1seq))
sdD_fx <- matrix(out$StanDev,nrow=length(x1seq))

zizoucols <- wes_palette("Zissou1",19,"continuous")
filled.contour(x1seq,x2seq,ED_fx,col= zizoucols,xlab=expression(O[i:j]),ylab=expression(I[i:j]),
               plot.axes = {axis(1);axis(2);points(x_j,pch=16)},lev=seq(min(ED_fx),max(ED_fx),len=20))

#true function plot
real_f <- matrix(500/(x[,1]-x[,2]),nrow=length(x1seq))
filled.contour(x1seq,x2seq,real_f,nlevels=40,color.palette=terrain.colors,xlab="x1",ylab="x2",
               plot.axes = {axis(1);axis(2);points(x_j,pch=16)})


#emulator sd
zizoucols2 <- wes_palette("Zissou1",27,"continuous")
filled.contour(x1seq,x2seq,sdD_fx,col=zizoucols2,xlab=expression(O[i:j]),ylab=expression(I[i:j]),
               plot.axes = {axis(1);axis(2);points(x_j,pch=16)})

#Define an observed liquidity ratio z, and the observation error sigma_e and model discrepancy sigma_epsilon
z <- 1.6				
sigma_e <- 0.04			#5% of |z| 
sigma_epsilon <- 0			#model discrepancy=0

#calculate the implausibility
Imp <- sqrt(  (ED_fx - z)^2 / ( sdD_fx^2 + sigma_epsilon^2 + sigma_e^2)  )

#plot the 2d implausibility (red high: input points are bad) 
levs <- c(seq(0,5,0.5),max(Imp))						#levels to colour the implausibility plot
zizoucols3 <- wes_palette("Zissou1",11,"continuous")

#plot the implausibility
filled.contour(x1seq,x2seq,Imp,col=zizoucols3,lev=levs,xlab=expression(O[i:j]),ylab=expression(I[i:j]),
               plot.axes = {axis(1);axis(2);points(x_j,pch=16); contour(x1seq, x2seq, Imp, levels = 3,lty = 1,lwd = 3, col = "white",labels = '3', add = TRUE)})

#run a wave 2
#select a set of 25 wave 2 runs "wave2" by choosing at random from the non-implausible points in x:
#LHS is smarter method to sample here
wave2 <- x[Imp<3,][sample(sum(Imp<3),25),]

x_j_wave1_2 <- rbind(x_j,wave2)
fx_j_wave2 <- 500/(wave2[,1]-wave2[,2])	
D_wave1_2 <- c(fx_j,fx_j_wave2)
out <- emulator_2d(x_j=x_j_wave1_2,x=x,D=D_wave1_2,beta_0=0,sig2=30,theta=c(55,55))
ED_fx <- matrix(out$Expect,nrow=length(x1seq))
sdD_fx <- matrix(out$StanDev,nrow=length(x1seq))

#We can plot the new emulator expectation over the 2d space:
filled.contour(x1seq,x2seq,ED_fx,col= zizoucols,xlab=expression(O[i:j]),ylab=expression(I[i:j]),
               plot.axes = {axis(1);axis(2);
                 points(x_j,pch=16)},lev=seq(min(ED_fx),max(ED_fx),len=20))


#calculate and plot the wave 2 implausibility
Imp_w2 <- sqrt(  (ED_fx - z)^2 / ( sdD_fx^2 + sigma_epsilon^2 + sigma_e^2)  )
filled.contour(x1seq,x2seq,Imp_w2,col=zizoucols2,lev=levs,xlab=expression(O[i:j]),ylab=expression(I[i:j]),
               plot.axes = {axis(1);axis(2);points(x_j,pch=16); 
                 points(x_j,pch=16);points(wave2,pch=16,col="yellow",cex=1);
                 contour(x1seq, x2seq, Imp_w2, levels = 3,lty = 1,lwd = 3, col = "white",labels = '3', add = TRUE)})

#proceed the same for further waves until non-implausible space is sufficiently reduced