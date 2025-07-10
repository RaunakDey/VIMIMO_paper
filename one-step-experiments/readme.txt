script_forward_problem -- solves the one step trajectory with SEIV 
15 minutes absoprtion and later no absoprtion model with Euler method

script_forward_problems_v2 -- same as before now with ode45, here initial conditions are
1e8 H and 1e7 virus and it just continues after dilution

all exposed theory -- use function one_step_all_exposed and v3
E1 has full storage of exposed host, rest is 0.

v4 -- just check proper adsorption rate.


simulator folder has one_step_simulate function to simulate, and before and after dilution function 
to simulate these conditions.

The time series of free and total phages are provided. Each time series has 3 columns as the experiments are performed in triplicates. 
Job is to quantify latent period and burst sizes



units: 'hrs'	'hrs'	'PFU/ml'	'PFU/ml'
variables: time	time	free phages		total phages