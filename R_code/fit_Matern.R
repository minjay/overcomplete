rm(list=ls())

library(RandomFields)
library(convoSPAT)
library(fields)
library(LANG)
library(INLA)
library(R.matlab)

RFoptions(modus_operandi="sloppy")

## read data
## change the directory if necessary
platform <- Sys.info()['sysname']
if (platform=='Darwin'){
  dat <- readMat('/Users/minjay/Documents/MATLAB/Needlets/overcomplete/real/data_sub.mat')
}else if (platform=='Linux'){
  dat <- readMat('/home/minjay/div_curl/MixedMatern/tests/real/wind.mat')
}

## re-arrange the data
samples <- t(dat$pot.samples)
n.obs <- nrow(samples)
theta <- dat$theta.samples
phi <- dat$phi.samples
theta <- dat$theta.samples*4

x <- sin(theta)*cos(phi)
y <- sin(theta)*sin(phi)  
z <- cos(theta)

## get distance matrix
## chordal distance
loc <- cbind(x, y, z)
Dist.mat <- as.vector(dist(loc))

#################################
## model definition            ##
#################################
mesh2 = inla.mesh.2d(loc = loc, max.edge = 0.1)
spde = inla.spde.create(mesh2, model = "matern")
