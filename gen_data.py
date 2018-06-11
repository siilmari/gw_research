from shared_stuff import *

# Let's define our model parameters A and alpha:
mu = [0.002, 0.00311]

parameters = [mu]*n_data

# Here, we're going to generate n_data datafiles into the folder 'data'.
for i in range(n_data):
	# First, we will calculate the expected EMRI rate for each bin, given the model parameters.
	# Bins will be arranged in a list in the following way: the first z_nbins elements will have
	# lnM=lnM_lower with z starting from z_lower and increasing all the way to z_upper-dz.
	# The next z_nbins elements will have lnM=lnM_lower+dlnM with z increasing, etc.
	rates_temp = list()
	for k in range(lnM_nbins):
		for j in range(z_nbins):
			rates_temp.append(rate(parameters[i], lnM_lower+k*dlnM, z_lower+j*dz))
	rates = np.multiply(rates_temp, volume_obslife_factors)

	# Now, we generate a vector of sample EMRI data using the rate vector:
	samples = np.random.poisson(rates)
	
	# Write the data into files:
	datafile = open("data/data%d.txt" % i, 'w+')
	datafile.write("The first ndim lines (after this one) contain the true parameters (A, alpha). (p019 is fixed to be "+str(p019)+".) The rest contain the rate of each bin and the corresponding sampled value (on the same line). Total number of events: "+str(sum(samples))+".\n")
	for k in range(ndim):
		datafile.write(str(parameters[i][k])+"\n") # write the parameters used to generate the data into the file
	for j in range(n_bins):
		datafile.write(str(rates[j])+" "+str(samples[j])+"\n")
	datafile.close()
