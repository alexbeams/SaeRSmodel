rm(list=ls()) #clears the working directory to start from a clean slate each time

# This code solves the SaeRS model in the manuscript using the invariant manifold.
# It simulates the nondimensional system

require(deSolve)	#load in the ODE solver package



###################
# external signal #
###################

#Assume that the external signal over time can be described by a piecewise linear function:

sgnlMax <- 1.4		#maximum for external signal, i.e. the degradation rate of P (modeling sequesestration of P as linear degradation/removal for simplicity)
tstart <- 2000000	#time when the ramping up commences
deltaT <- 1280000	#how long the signal requires to ramp to its maximum

#define the piecewise linear function describing the external signal
getSignal <- function(t) {
	if(t < tstart){0}else{
		if(t < tstart+deltaT){sgnlMax/deltaT*(t-tstart)}else{
			if(t < tstart+2*deltaT){sgnlMax - sgnlMax/deltaT*(t-tstart-deltaT)}else{0}
		}
	}
}
	

##############
# Parameters #		#description of parameter, [units of parameter]
##############


k1 <- 30   		#maximum rate at which P de-phosphorylates S* to S, [1/minutes]
k2 <- 80		#maximum rate at which Q de-phosphorylates S* to S, [1/minutes]
k3 <- 0.0085		#maximum rate at which S de-phosphorylates R* to R, [1/minutes]
k4 <- 0.0094		#maximum rate at which S* phosphorylates R to R*,   [1/minutes]
	
KP <- 10		#half-saturation constant for the enzyme P de-phosphorylating S, [concentration]
KQ <- 10		#half-saturation constant for the enzyme Q de-phosphorylating S, [concentration]
KRstar <- 10		#half-saturation constant for the enzyme S de-phosphorylating R*, [concentration]
KR <- 10		#half-saturation constant for the enyme S* phosphorylating R, [concentration]

kp3 <- 10		#baseline P3 promoter activity, [concentration/minute]
kp1 <- 42*kp3		#activated P1 promoter activity, [concentration/minute]

ks <- 90		#baseline rate of intrinsic phosphorylation of S, [1/minutes]
ksstar <- 5		#baseline rate of intrinsic de-phosphorylation of S*, [1/minutes]

beta <- .01		#mass-action constant for R* binding to free promoter site, [1/concentration/minutes]
k <- .8		#unbinding rate for R* molecules bound to promoter site, [1/minutes]

deltaR <- 1/9		#basal degradation/removal rate of R,R* molecules, [1/minutes]
deltaS <- 1/90		#basal degradation/removal rate of S,S* molecules, [1/minutes]
deltaP <- 1/90		#basal degradation/removal rate of P molecules, [1/minutes] 
deltaQ <- 1/40		#basal degradation/removal rate of Q molecules, [1/minutes]
deltaRstar <- deltaR

epsilon1 <- 1/90	#degradation rate for R* molecules when 1 promoter site is occupied, [1/minutes]
epsilon2 <- 1/10	#degradation rate for R* molecules when 2 promoter sites are occupied, [1/minutes]
epsilon3 <- 1/10	#degradation rate for R* molecules when 3 promoter sites are occupied, [1/minutes]
epsilon4 <- 1/10	#degradation rate for R* molecules when all 4 promoter sites are occupied, [1/minutes]
epsilon  <- epsilon1	#unequal values of these rates might result from stabilization/dimerization of R* molecules when bound to the promoter

xi	<- kp3/KP/2 #baseline production of SaeQ and SaeP (modeled deterministically, but ought to be stochastic)

#############
# timesteps #
#############

tEnd <- tstart+3*deltaT		#run the system with a slowly varying signal over a long period of time so we can ignore transient effects/visualize equilibrium behavior
DT <- 100			#temporal resolution for solver; can be rather coarse


##########################
# differential equations #
##########################

F <- function(Time, state, Pars) {
	with(as.list(c(state, Pars)), {


		chi4	<- z^4

		ds  	<- kp3/KP + kp1/KP*chi4 - ks*s - deltaS*s + ksstar*sigma + k1*sigma*p/(1+sigma) + k2*sigma*q/(KQ/KP+sigma) + k4*sigma*r/(1 + r) - k3*s*rho/(1+rho) 
		dsigma  <- -k1*sigma/(1 + sigma)*p-k2*sigma*q/(KQ/KP+sigma) + ks*s -(ksstar+deltaS)*sigma - k4*sigma*r/(1 + r)+k3*s*rho/(1+rho) 
		dr  	<- kp3/KR + kp1/KR*chi4 - deltaR*r + k3*s*KP/KR*rho/(1+rho) - k4*sigma*KP/KR*r/(1+r)
		drho  <-  -beta*rho*4*(1-z) + k*4*z/KRstar + k4*sigma*KP/KRstar*r/(1+r) - k3*s*KP/KRstar*rho/(1+rho) - deltaRstar*rho
		dp 	<- xi + kp1/KP*chi4-deltaP*p-p*getSignal(Time)
		dq	<- xi + kp1/KP*chi4-deltaQ*q
		dz	<- beta*rho*KRstar*(1-z) - (k+epsilon)*z	

	return(list(c(ds,dsigma,dr,drho,dp,dq,dz)))
	})
}



#########################
# set up and run solver #
#########################

pars <- c()				#defined parameters globally (sue me); feed in an empty "parameter" vector so the solver doesn't throw an error

Times <- seq(0, tEnd, by = DT)		#define timesteps

y0 <- c(s=kp3/deltaS/KP,sigma=0,r=kp3/deltaR/KR,rho=0.1/KRstar,p=kp3/deltaP/KP,q=kp3/deltaQ/KP,z=0)	#initial condition

out<- as.data.frame(ode(y0, Times, F, pars)) 		#solve the system! place solution in a dataframe called "out" (short for output)

out$S = out$s * KP					#undo the scaling; ensure we get the same solution as solving the original system
out$Sstar = out$sigma*KP
out$R = out$r*KR
out$Rstar = out$rho*KRstar
out$P = out$p * KP
out$Q = out$q * KP

out$sgnl <- sapply(out$time, getSignal )		#external signal is a time-dependent function fed into the ODE; obtain that signal for plotting

rbPal <- colorRampPalette(c('red','blue'))		#load a color palette for plotting bifurcation diagrams, if so desired

out$Col <- rbPal(1000)[as.numeric(cut(out$time,breaks = 1000))]		#this adds a column of color values based on the values

# thin out the times a little bit to make plotting easier
ind = seq(1,dim(out)[1],length=200)
ind = round(ind)
out=out[ind,]
out$time=out$time/1e6

#################
# plot solution #
#################


#setEPS()
#postscript("ramp_experiment.eps")

par(mfrow=c(3,2))		#make 4 plots in a 2x2 array

cxax=1.1
cxlb = 1.3

plot(sgnl~time,out,type='l',ylab=bquote(theta),xlab='time',lwd=2,cex.axis=cxax,cex.lab=cxlb)

plot(S~time,out,type='l',ylim=c(0,max(out$S,out$Sstar)),
	ylab='S,S*',lwd=2,cex.axis=cxax,cex.lab=cxlb)
lines(Sstar~time,out,lty='dashed',lwd=2)

plot(R~time,out,type='l',ylim=c(0,max(out$R,out$Rstar)),
	ylab='R,R*',lwd=2,cex.axis=cxax,cex.lab=cxlb)
lines(Rstar~time,out,lty='dashed',lwd=2)

plot(P~time,out,type='l',ylim=c(0,max(out$P)),ylab='P',lwd=2,cex.axis=cxax,cex.lab=cxlb)

plot(Q~time,out,type='l',ylim=c(0,max(out$Q)),ylab='Q',lwd=2,cex.axis=cxax,cex.lab=cxlb)

plot(z^4~time,out,type='l',ylim=c(0,1 ),ylab=bquote(chi[4]),lwd=2,cex.axis=cxax,cex.lab=cxlb)
	
#dev.off()

#make a bifurcation curve with direction arrows for time
#x11()
#par(mfrow=c(1,1))
#plot(Rstar~sgnl,out,type='l',xlab=bquote(theta),ylab=bquote(R^'*'))

inds = round( c(1/4 * dim(out)[1], 1/2 * dim(out)[1] ))
spc = 30

#for(arrow in 1:length(inds)){
	
#arrows( x0=out$sgnl[inds[arrow]-spc], y0 = out$Rstar[inds[arrow]-spc],
#		x1 = out$sgnl[inds[arrow]+spc], y1=out$Rstar[inds[arrow]+spc]  )

#}

