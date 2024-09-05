function dydt = one_step_eqn_after_dilution_infectionvary_absorbs(t,y,theta,NE)
%ONE_STEP_EQN Summary of this function goes here
%   Detailed explanation: Differential equations for SEIV (1,1,NE)

S = y(1);
E_mat(1:NE) = y(2:NE+1);
I = y(NE+2);
V_infec = y(NE+3);
V_noninfec = y(NE+4);

r = theta(1);
phi = theta(2);
tau = theta(3);
beta = theta(4);

prob_effective_infection = theta(6);
lysis_reset_coeff = theta(7);


etaeff = ((NE+1)/tau);

dotS = r*S - phi*V_infec*S*prob_effective_infection;
dotE1 = phi*S*V_infec*prob_effective_infection - etaeff * E_mat(1);
dotE_mat = etaeff.*E_mat(1:end-1) - etaeff.*E_mat(2:end);
dotI = etaeff * (E_mat(end) - I) ;
dotV_infec = (beta*prob_effective_infection) * etaeff * I - V_infec*phi*(S+I+sum(E_mat));
dotV_noninfec =  beta*(1-prob_effective_infection) * etaeff * I -  V_noninfec*phi*(S+I+sum(E_mat));

dydt = [dotS,dotE1,dotE_mat,dotI,dotV_infec,dotV_noninfec]';

end

