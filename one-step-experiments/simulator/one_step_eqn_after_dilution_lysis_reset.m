function dydt = one_step_eqn_after_dilution_lysis_reset(t,y,theta,NE)
%ONE_STEP_EQN Summary of this function goes here
%   Detailed explanation: Differential equations for SEIV (1,1,NE)

S = y(1);
E_mat(1:NE) = y(2:NE+1);
I = y(NE+2);
V_infec = y(NE+3);


r = theta(1);
phi = theta(2);
tau = theta(3);
beta = theta(4);

%prob_effective_infection = theta(6); %not used, everything is infectious
lysis_reset_coeff = theta(7);


etaeff = ((NE+1)/tau);

% all stages of reset is equal model
q(1:NE+1) = lysis_reset_coeff;


% % more reset at later stages model
% for i=1:NE+1
%     q(i) = lysis_reset_coeff*(i-1)/NE;
% end


%q(1) is not used 

dotS = r*S - phi*V_infec*S;
dotE1 = phi*S*V_infec- etaeff * E_mat(1) + sum(phi*V_infec*(q(2:NE).*E_mat(2:end)) ) + q(NE)*phi*I*V_infec;
dotE_mat = etaeff.*E_mat(1:end-1) - etaeff.*E_mat(2:end) - phi*V_infec*(q(2:NE).*E_mat(2:end));
dotI = etaeff * (E_mat(end) - I) - q(NE)*phi*I*V_infec;
dotV_infec = beta * etaeff * I - V_infec*phi*(S+I+sum(E_mat));

dydt = [dotS,dotE1,dotE_mat,dotI,dotV_infec]';

end


