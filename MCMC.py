import matplotlib.pyplot as pl
import emcee
import corner

from shared_stuff import *

# -----
# Here we will define the prior distribution, the likelihood the posterior distribution:
# -----

# The log prior. Outputs negative infinity (viewed as the "logarithm of zero")
# if the parameters lie outside the range of the prior (defined in the 'if' clause).
# The MCMC chain will then not venture into that area.
# --Input--
# p - vector of parameters
# --Output--
# The value of the prior at the specified parameter values
def lnprior(p):
	if p[0] > 0 and p[0] < 0.05 and p[1] > -1 and p[1] < 1:
	        return -m.log(p[0])
        else:
                return -np.inf

# The likelihood. Outputs negative infinity (viewed as the "logarithm of zero")
# if the parameters lie outside the range of the prior (defined in the 'if' clause).
# --Input--
# means - vector of intrinsic EMRI rates
# data - vector of values assumed to be from a Poisson distribution with means encoded in the vector 'means'
# --Output--
# The value of the likelihood at the specified parameter values
def lnlike(means, data):
	summa = 0.0
	for i in range(n_bins):
		if means[i] > 0:
			summa += -means[i] + data[i]*m.log(means[i])
	return summa

# The posterior distribution, which combines the prior and the likelihood. As before, it outputs negative
# infinity if the parameters lie outside the range of the prior.
# --Input--
# params - vector of parameters
# data - vector of values assumed to be from a Poisson distribution (the parameters found in the vector
#       'params' are used to calculate the means (i.e. the intrinsic EMRI rates) of the assumed distribution). 
# --Output--
# The value of the posterior at the specified parameter values
def lnprob(params, data):
	lnp = lnprior(params)
        if lnp == -np.inf:
                return -np.inf
        else:
		# First, we will use the parameter vector to calculate the corresponding intrinsic EMRI rates,
		# which we will assume to be the means of the underlying Poisson distribution:
		means_temp = list()
		for i in range(lnM_nbins):
			for j in range(z_nbins):
				means_temp.append(rate(params, lnM_lower+i*dlnM, z_lower+j*dz))
		means = np.multiply(means_temp, volume_obslife_factors)
		# Now, we use that to calculate the likelihood of the given parameter set:
                return lnp + lnlike(means, data)


# -----
# Now, we will move on to MCMC. First, we define some key variables:
# -----

# Set the number of "walkers" (different Markov Chains that run in parallel)
nwalkers = 100
# Number of runs (i.e. number of timesteps):
n_runs = 600
# The length of the burn-in period:
burn_in_len = 400

# Set a very crude initial guess for the parameters (used in generating the initial position of the walkers)
guess = [0.5, 0]

# -----
# Now, we will run MCMC, looping over all the datasets that need analysing.
# -----

for i in range(n_data):
	# First, we read in the data:
	print("Starting run %d" % i)
	datafile = open("data/data%d.txt" % i, 'r')
	datafile.readline()
	params = list()
	for j in range(ndim):
		params.append(float(datafile.readline().strip()))
	rates = list()
	samples = list()
	for k in range(n_bins):
		rate_str, sample_str = datafile.readline().split()
		rates.append(float(rate_str))
		samples.append(float(sample_str))
	datafile.close()

	# Let's initialize the walkers around our crude guess
	p0 = np.zeros([nwalkers, ndim])
	for j in range(len(p0)):
		p0[j] = np.random.multivariate_normal(guess, [[0.1,0], [0,0.2]])
		while lnprior(p0[j]) == -np.inf:
			# This step is taken so that the walkers wouldn't be initialised outside the "uniform prior zone" where the posterior is zero
			p0[j] = np.random.multivariate_normal(guess, [[0.1,0], [0,0.2]])

	# Now, let's create an EnsembleSampler object and sample with it (see the documentation of emcee for more details)
	sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args = [samples])
	print "Starting MCMC"
	pos, prob, state = sampler.run_mcmc(p0, n_runs)
	print "Finished MCMC"

	# Let's omit the first burn_in_len samples and make a flat Numpy array out of the rest of them:
	post_samples = sampler.chain[:, burn_in_len:, :].reshape((-1, ndim))

	# Now, let's find the element in the chain with the highest associated posterior probability:
	max_index = sampler.flatlnprobability.tolist().index(max(sampler.flatlnprobability.tolist()))
	max_lnprob = sampler.flatchain[max_index]
	print "The element in the chain with the highest associated posterior probability value: "+str(max_lnprob)

	# Now, let's find median estimates and construct the credible intervals for each dimension individually.
	# Let's stick them in one big Numpy array called CI (in each row we have the lower end of the CI, the
	# median and the higher end of the CI of each dimension):
	CI = np.empty([ndim, 3])
	for j in range(ndim):
		CI[j] = np.percentile(post_samples[:,j], [5, 50, 95])
	
	print("The median estimate for A is "+str(CI[0,1])+" and the 90% symmetric credible interval: ("+str(CI[0,0])+", "+str(CI[0,2])+")")
	print("The median estimate for alpha is "+str(CI[1,1])+" and the 90% symmetric credible interval: ("+str(CI[1,0])+", "+str(CI[1,2])+")")
	
	acceptance_frac = np.mean(sampler.acceptance_fraction)
	print "Mean acceptance fraction: "+str(acceptance_frac)

	# Let's write the results into a file:
	datafile = open("results/%d_estimate.txt" % i, 'w+')
	datafile.write("The next ndim lines contain the true parameters (A, alpha). (p019 is fixed to be "+str(p019)+")\n")
	for k in range(ndim):
		datafile.write(str(params[k])+"\n") # write the true parameters into the file
	datafile.write("The next ndim lines contain the parameters with the highest (combined) associated posterior probability (i.e. the mode):\n")
	for j in range(ndim):
		datafile.write(str(max_lnprob[j])+"\n") # write the estimated parameters into the file
	datafile.write("The next ndim lines contain the median estimates in each dimension:\n")
	for j in range(ndim):
		datafile.write(str(CI[j, 1])+"\n")
	datafile.write("The next ndim lines contain the 90% symmetric credible intervals in each dimension (each line contains the interval in one dimension: the lower and higher ends separated by a comma):\n")
	for j in range(ndim):
		datafile.write(str(CI[j, 0])+", "+str(CI[j, 2])+"\n")
	datafile.write("\n")
	datafile.write("Number of walkers: "+str(nwalkers)+"\n")
	datafile.write("Number of runs: "+str(n_runs)+"\n")
	datafile.write("Burn-in length: "+str(burn_in_len)+"\n")
	datafile.write("Mean acceptance fraction: "+str(acceptance_frac)+"\n")
	datafile.close()

	# Let's generate a corner plot. Note that some of the first iterations are not shown on the histograms
	# (we're essentially discarding the burn-in, making it pretty long to play it safe)
	fig_corner = corner.corner(post_samples, labels = [r"$A$", r"$\alpha$"], quantiles = [0.05, 0.5, 0.95], truths = params)
	fig_corner.savefig("results/%d_triangle_qs.png" % i)
	pl.close(fig_corner)

	# and a version without quantiles shown on the 1D histograms:
	fig_corner2 = corner.corner(post_samples, labels = [r"$A$", r"$\alpha$"], truths = params)
	fig_corner2.savefig("results/%d_triangle.png" % i)
	pl.close(fig_corner2)

	# Let's generate the trace plots
	trace_A = pl.figure()
	for j in range(nwalkers):
		pl.plot(sampler.chain[j, :, 0])
	pl.xlabel('Step #')
	pl.ylabel('Value')
	pl.title("Trace plot for the parameter A")
	trace_A.savefig("results/%d_traceplot_A.png" % i)
	pl.close(trace_A)

	trace_alpha = pl.figure()
	for j in range(nwalkers):
		pl.plot(sampler.chain[j, :, 1])
	pl.xlabel('Step #')
	pl.ylabel('Value')
	pl.title("Trace plot for the parameter alpha")
	trace_alpha.savefig("results/%d_traceplot_alpha.png" % i)
	pl.close(trace_alpha)

