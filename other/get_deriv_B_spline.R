library(fda)
library(R.matlab)

setwd('/Users/minjay/Documents/MATLAB/Needlets/overcomplete/real/prelim/')
data = readMat("theta_phi_R.mat")
x = data$theta[1, ]*4
breaks = c(0, 40/180, 80/180, 1)*pi

bS = bsplineS(x, breaks, nderiv=1)
writeMat('deriv_B_spline.mat', )