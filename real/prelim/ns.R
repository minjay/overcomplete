# this script obtains the natural cubic B-spline matrix evaluated at theta's
rm(list=ls())

require(splines)

library(R.matlab)

setwd('/Users/minjay/Documents/MATLAB/Needlets/overcomplete/real/prelim/')
data = readMat("theta_phi_R.mat")

xs = data$theta[1, ]*4
# intercept is not needed since it will be added explicitly
# out of 120/180*pi, the spline becomes linear
b_mat = ns(x=xs, knots=c(40/180*pi, 80/180*pi), intercept=FALSE, Boundary.knots=c(0, 120/180*pi))

matplot(xs, b_mat, type='l')

writeMat('ns.mat', b_mat=b_mat)
