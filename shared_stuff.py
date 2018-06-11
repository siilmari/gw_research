# This file contains the key variables and functions used in all other code files.
# All the variables, functions, etc. will be imported into all the other code files
# at the beginning.

import numpy as np
import math as m
from astropy.cosmology import Planck15

# Number of parameters to be estimated (A and alpha in this case):
ndim = 2

# Number of datasets to be generated and/or analysed:
n_data = 15

# Fix the parameter p_1:
p019 = -0.15

# Mass will be expressed in terms of 3 * 10^6 solar masses - define the scaling factor:
mass_scale_f = 3e6

# -----
# Now we will define the parameter ranges, bin widths, etc. of the log mass (lnM) and
# redshift (z) for the binned analysis:
# -----

lnM_lower = m.log(1e4/mass_scale_f) # lower end of the parameter range (approx 9.3)
lnM_upper = m.log(1e7/mass_scale_f) # upper end of the parameter range (approx 16.1)
lnM_nbins = 20 # number of bins
dlnM = (lnM_upper - lnM_lower)/lnM_nbins # bin width

# Similarly for z:
z_lower = 0.05
z_upper = 1.2
z_nbins = 20
dz = (z_upper - z_lower)/z_nbins

# Total number of bins:
n_bins = lnM_nbins*z_nbins


# -----
# Here we will define some functions and variables for defining the likelihood for our model:
# -----

# The main EMRI rate function, which calculates the intrinsic EMRI rate of a given bin.
# --Input--
# muu - vector of parameters
# lnM, z - lower ends of a particular bin
# --Output--
# The expected EMRI rate for the bin
def rate(muu, lnM, z):
	return r0(muu, lnM+dlnM/2) * muu[0] * m.exp(muu[1] * (lnM+dlnM/2)) * dlnM * dz
	# note that we add half of the bin width to lnM to calculate the rate for the midpoint of the bin
	# we will multiply with the comoving volume and observable lifetime separately (see volume_obslife_factors below)

# The R_0 factor for the EMRI rate function.
# --Input--
# muuu - vector of parameters
# lnM - log mass
# --Output--
# R_0
def r0(muuu, lnM):
	return 400*(3**p019)*m.exp(p019*lnM)
	# here, we should divide by 10e9 to get the rate in y^(-1) instead of Gyr^(-1);
	# however, we divide the differential comoving volume by that instead (and get it in Gpc instead of Mpc),
	# so it cancels out!

# The observable lifetime function - the result will be multiplied with the EMRI rate vector separately each time
# (see volume_obslife_factors below)
# --Input--
# lnM, z - log mass, redshift
# --Output--
# The observable lifetime (float)
def ObsLifeSpin4(lnM, z):
	zlist = [0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1,1.05,1.1,1.15,1.2,1.25,1.3,1.35,1.4,1.45,1.5,1.55,1.6,1.65,1.7,1.75,1.8,1.85,1.9,1.95,2,2.05,2.1,2.15,2.2,2.25,2.3,2.35,2.4,2.45,2.5,2.55,2.6,2.65,2.7,2.75,2.8,2.85,2.9,2.95,3]
	Mlist = [10000,12996.6,16891.1,21952.6,28530.8,37080.3,48191.6,62632.6,81400.9,105793,137495,178696,232243,301837,392284,509835,662610,861165,1.11922e+06,1.4546e+06,1.89048e+06,2.45698e+06,3.19323e+06,4.1501e+06,5.3937e+06,7.00996e+06,9.11054e+06,1.18406e+07,1.53887e+07,2e+07]
	obslife = [[23.8889,12.5758,6.66667,3.68056,1.93333,0.769231,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[22.1429,12.8788,7.3913,4.51389,2.73333,1.53846,0.740741,0.0595238,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[20.3968,12.5758,7.97101,5.13889,3.4,2.24359,1.48148,0.77381,0.287356,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[18.8095,11.9697,8.11594,5.625,3.86667,2.75641,1.97531,1.36905,0.862069,0.5,0.0537634,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[17.381,11.1364,8.04348,5.76389,4.2,3.14103,2.40741,1.84524,1.32184,1,0.645161,0.3125,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[16.1905,10.3788,7.6087,5.83333,4.46667,3.46154,2.71605,2.14286,1.72414,1.38889,1.07527,0.78125,0.505051,0.245098,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[15.1587,9.62121,7.17391,5.625,4.46667,3.58974,2.90123,2.38095,1.95402,1.61111,1.34409,1.09375,0.909091,0.686275,0.47619,0.277778,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[14.3651,8.93939,6.73913,5.34722,4.33333,3.52564,2.96296,2.5,2.12644,1.77778,1.55914,1.35417,1.16162,0.931373,0.809524,0.648148,0.45045,0.263158,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[13.8095,8.25758,6.23188,5,4.13333,3.46154,2.96296,2.5,2.18391,1.88889,1.66667,1.45833,1.31313,1.12745,1,0.87963,0.720721,0.614035,0.42735,0.25,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[13.3333,7.80303,5.7971,4.65278,3.86667,3.33333,2.83951,2.5,2.12644,1.88889,1.72043,1.51042,1.36364,1.22549,1.14286,1.01852,0.900901,0.789474,0.683761,0.541667,0.406504,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[12.8571,7.34848,5.36232,4.30556,3.6,3.14103,2.71605,2.38095,2.12644,1.88889,1.66667,1.51042,1.36364,1.27451,1.19048,1.06481,0.945946,0.877193,0.811966,0.666667,0.569106,0.436508,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[12.4603,7.04545,5.07246,4.02778,3.4,2.88462,2.53086,2.2619,2.01149,1.83333,1.6129,1.51042,1.36364,1.27451,1.14286,1.06481,0.990991,0.921053,0.854701,0.75,0.650407,0.515873,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[11.9841,6.74242,4.78261,3.75,3.13333,2.69231,2.40741,2.14286,1.89655,1.77778,1.6129,1.45833,1.31313,1.22549,1.19048,1.06481,0.990991,0.921053,0.811966,0.75,0.650407,0.0396825,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[11.5873,6.43939,4.49275,3.54167,2.93333,2.5,2.22222,1.96429,1.78161,1.66667,1.50538,1.35417,1.26263,1.17647,1.14286,1.01852,0.990991,0.877193,0.769231,0.583333,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[11.1905,6.06061,4.27536,3.26389,2.73333,2.30769,2.03704,1.84524,1.66667,1.55556,1.45161,1.30208,1.21212,1.12745,1.04762,0.972222,0.855856,0.526316,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[10.7937,5.60606,3.91304,3.05556,2.46667,2.11538,1.85185,1.66667,1.55172,1.44444,1.29032,1.19792,1.11111,1.02941,0.809524,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[10.2381,5.15152,3.55072,2.77778,2.26667,1.92308,1.7284,1.54762,1.43678,1.27778,1.12903,0.885417,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[9.60317,4.62121,3.18841,2.43056,2,1.73077,1.54321,1.36905,1.03448,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[8.88889,4.09091,2.68116,2.08333,1.66667,1.34615,0.246914,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[7.93651,3.40909,2.24638,1.59722,0.4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[6.8254,2.65152,1.30435,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[5.39683,1.21212,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[3.65079,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[1.34921,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]]

	NM=30
	Nz=60	
	M = mass_scale_f * m.exp(lnM)

	i=0
	if M < Mlist[0]:
		i = 1
		print "Warning - mass below minimal mass in function!"
	elif M > Mlist[NM-1]:
		i = NM - 1
		print "Warning - mass above maximal mass in function!"
	else:
		while Mlist[i] <= M and i < NM-1:
			i += 1
	j = 0
	if z < zlist[0]:
		j = 1
		print "Warning - redshift below minimal redshift in function!"
	elif z > zlist[Nz-1]:
		j = Nz - 1
		print "Warning - redshift above maximal redshift in function!"
	else:
		while zlist[j] <= z and j < Nz-1:
			j += 1

	life = obslife[i-1][j-1] * (Mlist[i] - M) * (zlist[j] - z)
	life += obslife[i][j-1] * (M-Mlist[i-1]) * (zlist[j] - z)
	life += obslife[i-1][j] * (Mlist[i]-M) * (z - zlist[j-1])
	life += obslife[i][j] * (M-Mlist[i-1]) * (z - zlist[j-1])
	life /= (Mlist[i] - Mlist[i-1]) * (zlist[j] - zlist[j-1])
	return life

# Generate a vector of values which will be multiplied with the corresponding EMRI rate in each bin.
# These factors represent the differential comoving volumes and observable lifetimes.
volume_obslife_factors = list()
for k in range(lnM_nbins):
	for j in range(z_nbins):
		lnM_eff = lnM_lower + k*dlnM + dlnM/2 # add half of the bin width to get the midpoint
		z_eff = z_lower + j*dz + dz/2
		volume_obslife_factors.append(ObsLifeSpin4(lnM_eff, z_eff) * Planck15.differential_comoving_volume(z_eff).value / 10e9 * 4 * m.pi)
