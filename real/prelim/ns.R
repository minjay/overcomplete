# this script obtains the natural cubic B-spline matrix evaluated at theta's
rm(list=ls())

require(splines)

library(R.matlab)
library(JM)

setwd('/Users/minjay/Documents/MATLAB/Needlets/overcomplete/real/prelim/')
data = readMat("theta_phi_R.mat")

xs = data$theta[1, ]*4
# intercept is not needed since it will be added explicitly
# out of 120/180*pi, the spline becomes linear
b_mat = ns(x=xs, knots=c(40/180*pi, 80/180*pi), intercept=FALSE, Boundary.knots=c(0, 120/180*pi))

b_mat_fine = ns(x=xs, knots=c(30/180*pi, 60/180*pi, 90/180*pi), intercept=FALSE, Boundary.knots=c(0, 120/180*pi))

# derivatives of basis functions at the same locations
b_mat_deriv = dns(x=xs, knots=c(40/180*pi, 80/180*pi), intercept=FALSE, Boundary.knots=c(0, 120/180*pi))

matplot(xs, b_mat, type='l')

matplot(xs, b_mat_deriv, type='l')

writeMat('ns.mat', b_mat=b_mat, b_mat_fine=b_mat_fine)

writeMat('ns_deriv.mat', b_mat_deriv=b_mat_deriv)
