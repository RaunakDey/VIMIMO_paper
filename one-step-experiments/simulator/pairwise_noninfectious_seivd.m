function dydt = pairwise_noninfectious_seivd(t,y,theta,NE,Dc)
% ONE_STEP_EQN Summary of this function goes here
% Detailed explanation: Differential equations for SEIV (1,1,NE)


S = y(1);
E_mat(1:NE) = y(2:NE+1);
I = y(NE+2);
V_I = y(NE+3);
V_NI = y(NE+4);
D = y(NE+5);

r = theta(1);
phi = theta(2);
tau = theta(3);
beta = theta(4);

etaeff = ((NE+1)/tau);

hill_factor = 1/(1+(D/Dc)^2);

dotS = r*S - phi*V_I*S*hill_factor ;
dotE1 = phi*S*V_I*hill_factor - etaeff * E_mat(1);
dotE_mat = etaeff.*E_mat(1:end-1) - etaeff.*E_mat(2:end);
dotI = etaeff * (E_mat(end) - I);
dotV_I = beta * etaeff * I - hill_factor*V_I*phi*(S+I+sum(E_mat));
dotV_NI = 0;
dotD =  etaeff * I;

dydt = [dotS,dotE1,dotE_mat,dotI,dotV_I,dotV_NI,dotD]';

end


