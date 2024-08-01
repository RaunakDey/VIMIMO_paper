<p align="center">
  <picture>
<img width="1417" alt="logo" src="https://github.com/RaunakDey/VIMIMO-Virus-Microbe-modeller/assets/39820997/2a441485-ccb7-445d-917c-bd76aa658daf">
  </picture>
</p>


<h1 align="center"><em>Scaling phage–bacteria dynamics from pairwise interactions to complex communities</em>: Code repository</h1>

*By Raunak Dey, 2024*


This is the repository for the paper -- if you want to take a more detailed look you may find our project [VIMIMO](https://github.com/RaunakDey/VIMIMO-Virus-Microbe-modeller/) useful. In this repository all the ODE models, the Bayesian scripts and figure generating scripts can be found. The Bayesian analysis is performed used DRAM-MCMC algorithm using the open-sourced mcmcstat package. 

## Requirements and installation

You will need MATLAB 2022b and higher to run this code. You need to install the mcmcstat package and add the directory to the matlab path which can be done with.
```matlab
git clone git@github.com:mjlaine/mcmcstat.git
addpath(genpath(./mcmcstat))
```
The data files and the mcmc chains come with the directory for the onesteps. For the community analysis you need to download the datafiles and put it in the proper directory. The results folder should be placed inside 
```
./community/
```


## Notes for the code reviewer:
We invite all reviewers to explore our models. However here is a suggested outline if you want to follow along with the paper.
- The onesteps are modeled using SEIV model, the ODE for that can be found in SI-1.1 in the paper and the code can be found with the [ODE models](./one-step-experiments/simulator/one_step_eqn_before_dilution.m) and [simulation function](./one-step-experiments/simulator/one_step_eqn_before_dilution.m) one_step_simulate_particular_points.m takes care of the the dilution.
- For the Bayesian analysis of each of the interactions we can look at the script_inverse_interactionsi.m files where $i=1,2,\dots,8$ for the 8 onestep growth curves. The one-step data is saved in [data_2024](./one-step-experiments/data_2024/) and the mcmc run results are stored in [result_replicates](./one-step-experiments/result_replicates/).
- For the community experiments, the data is stored [here](./community/data/) and the mcmcruns are being downloaded.
- We invite the reviewers to check the SEIV and SEIVD models in SI-1.3 and SI-1.4 respectively in the paper. The code for the same can be found in 

