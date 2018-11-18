# Single-host species preference model
# Preference depends on vector status (non-viruliferous versus viruliferous)
# Lauren Shoemaker

# Requires library deSolve to solve differential equations
# Load library 
library(deSolve)

# Model set-up
times <- 0:100            # Total amount of time-steps to run the model
T <- 100                  # Total number of plants (won't alter outcome)

starting_infection <- .1  # Percent hosts already infected
starting_vectors <- 50    # Starting vector population size

sP <- round((1-starting_infection)*T)  # Starting number fo healthy hosts
iP <-round(starting_infection*T)       # Starting number of infected hosts
sV <- starting_vectors                 # Starting number of healthy vectors
iV <- 0                                # Starting number of infected vectors
starting <- c(sP, iP, sV, iV)

# Model parameters
r_I <- .2245              # Vector growth on infected hosts
r_H <- .2245              # Vector growth on healthy hosts
K <- 5                    # Carrying capacity 
beta_VP <- .05            # Transmission vector to host
beta_PV <- .05            # Transmission host to vector

# ----------------------------------------------------------------------------------------------------
# Create system of differential equations that will be solved using lsoda
system.eqns=function(t,N,pars) {
  
  # total number of equations
  Ndot <- numeric(4); 
  
  # susceptible host 
  Ndot[1] = -beta_VP*N[4]*(I_H*N[1]/(I_H*N[1]+I_I*N[2]))
  
  # infected host 
  Ndot[2] = beta_VP*N[4]*(I_H*N[1]/(I_H*N[1]+I_I*N[2]))
  
  # susceptible vector
  Ndot[3] = r_H*((H_H*N[1]/(H_H*N[1]+H_I*N[2]))*N[3] + 
                   (I_H*N[1]/(I_H*N[1]+I_I*N[2]))*N[4])*(1-((H_H*N[1]/(H_H*N[1]+H_I*N[2]))*N[3] + 
                                                              (I_H*N[1]/(I_H*N[1]+I_I*N[2]))*N[4])/(K*N[1])) - 
    beta_PV*(1 - (H_H*N[1]/(H_H*N[1]+H_I*N[2])))*N[3]
  
  # infected vector
  Ndot[4] = r_I*((1 - (H_H*N[1]/(H_H*N[1]+H_I*N[2])))*N[3] + 
                   (1 - (I_H*N[1]/(I_H*N[1]+I_I*N[2])))*N[4])*(1-((1 - (H_H*N[1]/(H_H*N[1]+H_I*N[2])))*N[3] + 
                                                                    (1 - (I_H*N[1]/(I_H*N[1]+I_I*N[2])))*N[4])/(K*N[2])) + 
    beta_PV*(1 - (H_H*N[1]/(H_H*N[1]+H_I*N[2])))*N[3]
  
  return(list(Ndot));
  
}

# -------------------------------------------------------------------------------------------------
# Scenarios to Run

# 1 No preference
H_I <- .5 # Preference of healthy vectors for infected hosts
I_I <- .5 # Preference of infected vectors for infected hosts
H_H <- .5 # Preference of healthy vectors for healthy hosts
I_H <- .5 # Preference of infected vectors for healthy hosts

pars <- c(r_H, r_I, K, beta_VP, beta_PV, H_I, H_H, I_I, I_H, T) # Concatenate parameters

# Solve using lsoda
soln_no=ode(y=starting, times=times, func=system.eqns, parms=pars, hmax=0.1)

# 2 Full opposite preference
H_I <- .85 # Preference of healthy vectors for infected hosts
I_I <- .15 # Preference of infected vectors for infected hosts
H_H <- .15 # Preference of healthy vectors for healthy hosts
I_H <- .85 # Preference of infected vectors for healthy hosts

pars <- c(r_H, r_I, K, beta_VP, beta_PV, H_I, H_H, I_I, I_H, T) # Concatenate parameters

# Solve using lsoda
soln_full_opp=lsoda(y=starting, times=times, func=system.eqns, parms=pars, hmax=0.01)

# 3 Full same preference
H_I <- .15 # Preference of healthy vectors for infected hosts
I_I <- .85 # Preference of infected vectors for infected hosts
H_H <- .85 # Preference of healthy vectors for healthy hosts
I_H <- .15 # Preference of infected vectors for healthy hosts

pars <- c(r_H, r_I, K, beta_VP, beta_PV, H_I, H_H, I_I, I_H, T) # Concatenate parameters

# Solve using lsoda
soln_full=lsoda(y=starting, times=times, func=system.eqns, parms=pars, hmax=0.1)

# 4 Observed  preference
H_I <- .63 # Preference of healthy vectors for infected hosts
I_I <- .49 # Preference of infected vectors for infected hosts
H_H <- .37 # Preference of healthy vectors for healthy hosts
I_H <- .51 # Preference of infected vectors for healthy hosts

pars <- c(r_H, r_I, K, beta_VP, beta_PV, H_I, H_H, I_I, I_H, T) # Concatenate parameters

# Solve using lsoda
soln_observed=lsoda(y=starting, times=times, func=system.eqns, parms=pars, hmax=0.1)

# ----------------------------------------------------------------------------------------------------
# Example Plot

# Plotting colors, color-blind friendly
red <- rgb(213,94,0,max=255)
blue <- rgb(86,180,233,max=255)
green <- rgb(0,158,115,max=255)

# Plot results from model run
quartz(width=9, height=4)
par(mfrow=c(1,3))

# First plot host population results
plot(soln_no[,3], type="l", col="grey", lwd=3, ylim=c(0,100), ylab="Host Population (% Infected)", 
     xlab="Time", xaxt="n", yaxt="n")
lines(soln_full[,3], col=green, lwd=3, lty=2)
lines(soln_full_opp[,3], col=red, lwd=3, lty=2)
lines(soln_observed[,3], col="black", lwd=3, lty=4)
axis(side=1, at=c(0, 50, 100), tick=TRUE, labels=TRUE)
axis(side=1, at=c(0, 25, 50, 75, 100), tick=TRUE, labels=FALSE)
axis(side=2, at=c(0, 25, 50, 75, 100), tick=TRUE, labels=FALSE)
axis(side=2, at=c(0, 50, 100), tick=TRUE, labels=TRUE)

# Then plot vector population results
plot(100*soln_no[,5]/(soln_no[,4]+soln_no[,5]), type="l", col="grey", lwd=3, ylim=c(0,100), 
     ylab="Vector Population (% Infected)", xlab="Time", xaxt="n", yaxt="n")
lines(100*soln_full[,5]/(soln_full[,4]+soln_full[,5]), col=green, lwd=3, lty=2)
lines(100*soln_full_opp[,5]/(soln_full_opp[,4]+soln_full_opp[,5]), col=red, lwd=3, lty=2)
lines(100*soln_observed[,5]/(soln_observed[,4]+soln_observed[,5]), col="black", lwd=3, lty=4)
axis(side=1, at=c(0, 50, 100), tick=TRUE, labels=TRUE)
axis(side=1, at=c(0, 25, 50, 75, 100), tick=TRUE, labels=FALSE)
axis(side=2, at=c(0, 25, 50, 75, 100), tick=TRUE, labels=FALSE)
axis(side=2, at=c(0, 50, 100), tick=TRUE, labels=TRUE)

# Add legend
plot.new()
legend("center", c("No Preference ", "Opposite Preference","Same Preference", "Observed Preference"), 
       col=c("grey",red, green, "black"), 
       lwd=c(3,3,3,3), lty=c(1, 2, 2, 4), bty="n", seg.len = 3.5)
