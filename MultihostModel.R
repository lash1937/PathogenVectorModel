# Multihost species preference model
# Preference depends on vector status (non-viruliferous versus viruliferous) as in SingleHostModel.R
# And now also depends on species identity.
# This code is written for two species identities: annual versus perennial hosts.
# Lauren Shoemaker

# Requires library deSolve to solve differential equations
# Load library 
library(deSolve)

# Model set-up
times <- 0:100            # Total amount of time-steps to run the model
T <- 100                  # total number of plants
                          # (won't alter outcome if starting with a percent of hosts infected)
percent_annual <- .5      # Percent of the host community consisting 

starting_infection <- .1  # Percent hosts already infected
starting_vectors <- 50    # Starting vector population size


sA <- round(T*percent_annual*(1-starting_infection))        # Starting number of healthy annuals
iA <- round(T*percent_annual*starting_infection)            # Starting number of infected annuals
sP <- round(T*(1-percent_annual)*(1-starting_infection))    # Starting number of healthy perennials
iP <-round(T*(1-percent_annual)*starting_infection)         # Starting number of infected perennials
sV <- starting_vectors                                      # Starting number of healthy vectors
iV <- 0                                                     # Starting number of infected vectors
starting <- c(sA, iA, sP, iP, sV, iV)

# Model parameters
r_AI <- .2245       # Vector growth on infected annual
r_AH <- .2245       # Vector growth on healthy annual
r_PI <- .2245       # Vector growth on infected perennial
r_PH <- .2245       # Vector growth on healthy perennial
K_A <- 5            # Carrying capacity annual
K_P <- 5            # Carrying capacity perennial
beta_VA <- .05      # Transmission vector to annual
beta_VP <- .05      # Transmission vector to perennial
beta_AV <- .05      # Transmission annual to vector
beta_PV <- .05      # Transmission perennial to vector


# ----------------------------------------------------------------------------------------------------
# Create system of differential equations that will be solved using lsoda
system.eqns=function(t,N,pars) {
  
  # total number of equations
  Ndot <- numeric(6); 

  # susceptible annual 
  Ndot[1] = -beta_VA*N[6]*(i_AH*N[1]/(i_AH*N[1]+i_AI*N[2]+i_PH*N[3]+i_PI*N[4]))
  
  # infected annual 
  Ndot[2] = beta_VA*N[6]*(i_AH*N[1]/(i_AH*N[1]+i_AI*N[2]+i_PH*N[3]+i_PI*N[4]))
  
  # susceptible perennial 
  Ndot[3] = -beta_VP*N[6]*(i_PH*N[3]/(i_AH*N[1]+i_AI*N[2]+i_PH*N[3]+i_PI*N[4]))
  
  # infected perennial 
  Ndot[4] = beta_VP*N[6]*(i_PH*N[3]/(i_AH*N[1]+i_AI*N[2]+i_PH*N[3]+i_PI*N[4]))
  
  # susceptible vector
  Ndot[5] = r_AH*((h_AH*N[1]/(h_AH*N[1]+h_AI*N[2]+h_PH*N[3]+h_PI*N[4]))*
                    N[5]+(i_AH*N[1]/(i_AH*N[1]+i_AI*N[2]+i_PH*N[3]+i_PI*N[4]))*N[6])*
    (1-(((h_AH*N[1]/(h_AH*N[1]+h_AI*N[2]+h_PH*N[3]+h_PI*N[4]))*
           N[5]+(i_AH*N[1]/(i_AH*N[1]+i_AI*N[2]+i_PH*N[3]+i_PI*N[4]))*N[6])/(K_A*N[1]))) + 
    r_PH*((h_PH*N[3]/(h_AH*N[1]+h_AI*N[2]+h_PH*N[3]+h_PI*N[4]))*N[5]+
            (i_PH*N[3]/(i_AH*N[1]+i_AI*N[2]+i_PH*N[3]+i_PI*N[4]))*N[6])*
    (1-(((h_PH*N[3]/(h_AH*N[1]+h_AI*N[2]+h_PH*N[3]+h_PI*N[4]))*N[5]+
           (i_PH*N[3]/(i_AH*N[1]+i_AI*N[2]+i_PH*N[3]+i_PI*N[4]))*N[6])/(K_P*N[3]))) - 
    beta_AV*(h_AI*N[2]/(h_AH*N[1]+h_AI*N[2]+h_PH*N[3]+h_PI*N[4]))*N[5] - 
    beta_PV*(h_PI*N[4]/(h_AH*N[1]+h_AI*N[2]+h_PH*N[3]+h_PI*N[4]))*N[5]
  
  # infected vector
  Ndot[6] = r_AI*((h_AI*N[2]/(h_AH*N[1]+h_AI*N[2]+h_PH*N[3]+h_PI*N[4]))*
                    N[5]+(i_AI*N[2]/(i_AH*N[1]+i_AI*N[2]+i_PH*N[3]+i_PI*N[4]))*N[6])*
    (1-(((h_AI*N[2]/(h_AH*N[1]+h_AI*N[2]+h_PH*N[3]+h_PI*N[4]))*
           N[5]+(i_AI*N[2]/(i_AH*N[1]+i_AI*N[2]+i_PH*N[3]+i_PI*N[4]))*N[6])/(K_A*N[2]))) + 
    r_PI*((h_PI*N[4]/(h_AH*N[1]+h_AI*N[2]+h_PH*N[3]+h_PI*N[4]))*N[5]+
            (i_PI*N[4]/(i_AH*N[1]+i_AI*N[2]+i_PH*N[3]+i_PI*N[4]))*N[6])*
    (1-(((h_PI*N[4]/(h_AH*N[1]+h_AI*N[2]+h_PH*N[3]+h_PI*N[4]))*N[5]+
           (i_PI*N[4]/(i_AH*N[1]+i_AI*N[2]+i_PH*N[3]+i_PI*N[4]))*N[6])/(K_P*N[4]))) + 
    beta_AV*(h_AI*N[2]/(h_AH*N[1]+h_AI*N[2]+h_PH*N[3]+h_PI*N[4]))*N[5] + 
    beta_PV*(h_PI*N[4]/(h_AH*N[1]+h_AI*N[2]+h_PH*N[3]+h_PI*N[4]))*N[5]
  
  return(list(Ndot));
  
}

# -------------------------------------------------------------------------------------------------
# Scenarios to Run

# 1 No preference

H_I <- .5 # Preference of healthy aphids for infected plants
I_I <- .5 # Preference of infected aphids for infected plants
H_H <- .5 # Preference of healthy aphids for healthy plants
I_H <- .5 # Preference of infected aphids for healthy plants
H_P <- .5 # Preference of healthy aphids for perennial plants
H_A <- .5 # Preference of healthy aphids for annual plants
I_P <- .5 # Preference of infected aphids for perennial plants
I_A <- .5 # Preference of infected aphids for annual plants

# Composite probabilities
h_AI <- H_A*H_I
h_AH <- H_A*H_H
h_PI <- H_P*H_I
h_PH <- H_P*H_H
i_AI <- I_A*I_I
i_AH <- I_A*I_H
i_PI <- I_P*I_I
i_PH <- I_P*I_H

pars <- c(r_AH, r_AI, r_PH, r_PI, K_A, K_P, beta_VP, beta_PV, beta_AV, beta_VA, 
          h_AI, h_AH, h_PI, h_PH, i_AI, i_AH, i_PI, i_PH,  T)

# Solve system of differential equations
soln_no=ode(y=starting, times=times, func=system.eqns, parms=pars, hmax=0.1)

# 2 Specialist preference

H_I <- .5  # Preference of healthy aphids for infected plants
I_I <- .5  # Preference of infected aphids for infected plants
H_H <- .5  # Preference of healthy aphids for healthy plants
I_H <- .5  # Preference of infected aphids for healthy plants
H_P <- .85 # Preference of healthy aphids for perennial plants
H_A <- .15 # Preference of healthy aphids for annual plants
I_P <- .85 # Preference of infected aphids for perennial plants
I_A <- .15 # Preference of infected aphids for annual plants

# Composite probabilities
h_AI <- H_A*H_I
h_AH <- H_A*H_H
h_PI <- H_P*H_I
h_PH <- H_P*H_H
i_AI <- I_A*I_I
i_AH <- I_A*I_H
i_PI <- I_P*I_I
i_PH <- I_P*I_H

# Concatenate parameters
pars <- c(r_AH, r_AI, r_PH, r_PI, K_A, K_P, beta_VP, beta_PV, beta_AV, beta_VA, 
          h_AI, h_AH, h_PI, h_PH, i_AI, i_AH, i_PI, i_PH,  T)

# Solve system of differential equations
soln_s=ode(y=starting, times=times, func=system.eqns, parms=pars, hmax=0.1)

# 3 Observed preference for annuals versus perennials

H_I <- .5 # preference of healthy aphids for infected plants
I_I <- .5 # preference of infected aphids for infected plants
H_H <- .5 # preference of healthy aphids for healthy plants
I_H <- .5 # preference of infected aphids for healthy plants
H_P <- .43 # preference of healthy aphids for perennial plants
H_A <- .57 # preference of healthy aphids for annual plants
I_P <- .79 # preference of infected aphids for perennial plants
I_A <- .21 # preference of infected aphids for annual plants

# Composite probabilities
h_AI <- H_A*H_I
h_AH <- H_A*H_H
h_PI <- H_P*H_I
h_PH <- H_P*H_H
i_AI <- I_A*I_I
i_AH <- I_A*I_H
i_PI <- I_P*I_I
i_PH <- I_P*I_H

# Concatenate parameters
pars <- c(r_AH, r_AI, r_PH, r_PI, K_A, K_P, beta_VP, beta_PV, beta_AV, beta_VA, 
          h_AI, h_AH, h_PI, h_PH, i_AI, i_AH, i_PI, i_PH,  T)

# Solve system of differential equations
soln_obs=ode(y=starting, times=times, func=system.eqns, parms=pars, hmax=0.1)

# 4 Specialist to generalist switch

H_I <- .5 # preference of healthy aphids for infected plants
I_I <- .5 # preference of infected aphids for infected plants
H_H <- .5 # preference of healthy aphids for healthy plants
I_H <- .5 # preference of infected aphids for healthy plants
H_P <- .85 # preference of healthy aphids for perennial plants
H_A <- .15 # preference of healthy aphids for annual plants
I_P <- .5 # preference of infected aphids for perennial plants
I_A <- .5 # preference of infected aphids for annual plants

# Composite probabilities
h_AI <- H_A*H_I
h_AH <- H_A*H_H
h_PI <- H_P*H_I
h_PH <- H_P*H_H
i_AI <- I_A*I_I
i_AH <- I_A*I_H
i_PI <- I_P*I_I
i_PH <- I_P*I_H

# Concatenate parameters
pars <- c(r_AH, r_AI, r_PH, r_PI, K_A, K_P, beta_VP, beta_PV, beta_AV, beta_VA, 
          h_AI, h_AH, h_PI, h_PH, i_AI, i_AH, i_PI, i_PH,  T)

# Solve system of differential equations
soln_s_to_g=ode(y=starting, times=times, func=system.eqns, parms=pars, hmax=0.1)

# 5 Generalist to specialist switch

H_I <- .5 # preference of healthy aphids for infected plants
I_I <- .5 # preference of infected aphids for infected plants
H_H <- .5 # preference of healthy aphids for healthy plants
I_H <- .5 # preference of infected aphids for healthy plants
H_P <- .5 # preference of healthy aphids for perennial plants
H_A <- .5 # preference of healthy aphids for annual plants
I_P <- .85 # preference of infected aphids for perennial plants
I_A <- .15 # preference of infected aphids for annual plants

# Composite probabilities
h_AI <- H_A*H_I
h_AH <- H_A*H_H
h_PI <- H_P*H_I
h_PH <- H_P*H_H
i_AI <- I_A*I_I
i_AH <- I_A*I_H
i_PI <- I_P*I_I
i_PH <- I_P*I_H

# Concatenate parameters
pars <- c(r_AH, r_AI, r_PH, r_PI, K_A, K_P, beta_VP, beta_PV, beta_AV, beta_VA, 
          h_AI, h_AH, h_PI, h_PH, i_AI, i_AH, i_PI, i_PH,  T)

# Solve system of differential equations
soln_g_to_s=ode(y=starting, times=times, func=system.eqns, parms=pars, hmax=0.1)

# ----------------------------------------------------------------------------------------------------
# Example Plot for a community of 50/50 annuals and perennials
# Recreates figure 3c and 3d from the manuscript

# Plotting colors, color-blind friendly
red <- rgb(213,94,0,max=255)
green <- rgb(0,158,115,max=255)
blue <- rgb(86,180,233,max=255)

# Plot results from model run
quartz(width=9, height=4)
par(mfrow=c(1,3))

# First plot host population results
plot(soln_no[,3] + soln_no[,5], type="l", col="grey", lwd=3, ylim=c(0,100), 
     ylab="Host Population (% Infected)", xlab="Time", xaxt="n", yaxt="n")
axis(side=1, at=c(0, 25, 50, 75, 100), tick=TRUE, labels=TRUE)
axis(side=2, at=c(0, 25, 50, 75, 100), tick=TRUE, labels=FALSE)
axis(side=2, at=c(0, 50, 100), tick=TRUE, labels=TRUE)
lines(soln_s_to_g[,3] + soln_s_to_g[,5], col=green, lwd=3, lty=2)
lines(soln_g_to_s[,3] + soln_g_to_s[,5], col=red, lwd=3, lty=2)
lines(soln_s[,3] + soln_s[,5], col=blue, lwd=3, lty=1)
lines(soln_obs[,3] + soln_obs[,5], col="black", lwd=3, lty=4)

# Then plot vector population results
plot(100*soln_no[,7]/(soln_no[,6]+soln_no[,7]), type="l", col="grey", 
     lwd=3, ylim=c(0,100), ylab="Vector Population (% Infected)", xlab="Time", xaxt="n", yaxt="n")
axis(side=1, at=c(0, 25, 50, 75, 100), tick=TRUE, labels=TRUE)
axis(side=2, at=c(0, 25, 50, 75, 100), tick=TRUE, labels=FALSE)
lines(100*soln_s_to_g[,7]/(soln_s_to_g[,6]+soln_s_to_g[,7]), col=green, lwd=3, lty=2)
lines(100*soln_g_to_s[,7]/(soln_g_to_s[,6]+soln_g_to_s[,7]), col=red, lwd=3, lty=2)
lines(100*soln_s[,7]/(soln_s[,6]+soln_s[,7]), col=blue, lwd=3, lty=1)
lines(100*soln_obs[,7]/(soln_obs[,6]+soln_obs[,7]), col="black", lwd=3, lty=4)

plot.new()
legend("center", c("Generalist", "Specialist", 
                   "Generalist to Specialist Preference",
                   "Specialist to Generalist Preference", 
                   "Observed Preference"), 
       col=c("grey",blue, red, green, "black"), 
       lwd=c(3,3, 3,3,3), lty=c(1, 1, 2, 2, 4), bty="n", seg.len = 3.5)
