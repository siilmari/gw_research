----------------------
Introduction
----------------------

This is the code I used for analysing gravitational wave data for my Bachelor's Thesis. The aim of the project was to assess the potential of LISA, the space-based gravitational wave detector being built by the European Space Agency, for studying the evolution and mass distribution of black holes. Specifically, this project concentrated on extreme mass ratio inspirals (EMRIs), astronomical systems comprised of a stellar-mass black hole spiralling into a massive black hole. It was done by generating sample EMRI data and analysing it using Markov Chain Monte Carlo (MCMC). The data analysis was implemented in Python, using the package emcee.

For the project, many different versions of the code were used, tailored to each specific application (say, estimating different sets of variables). This repository contains one such version of the code.

For a more detailed explanation of the astrophysical background and data analysis methods, see my Bachelor's Thesis, which is included in the repository. The results and future work are also discussed there.

----------------------
Structure and usage
----------------------

* shared_stuff.py - the functions and variables that are used by all other code files
* gen_data.py - LISA data generation
* MCMC.py - data analysis

The code is entirely self-contained: one can generate sample data and analyse it without needing any additional files - try it! First, two folders need to be created in the directory: 'data' and 'results'. Then, run 'gen_data.py', which will save some data files into the folder 'data'. After that, by running 'MCMC.py', these datasets will be analysed. The resulting parameter estimates and plots will be saved into the folder 'results'. The Thesis will explain what they mean.