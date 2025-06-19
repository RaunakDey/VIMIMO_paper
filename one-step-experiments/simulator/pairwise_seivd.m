function dydt = pairwise_seivd(t,y,theta,NE,Dc)
%ONE_STEP_EQN Summary of this function goes here
%   Detailed explanation: Differential equations for SEIV (1,1,NE)


S = y(1);
E_mat(1:NE) = y(2:NE+1);
I = y(NE+2);
V = y(NE+3);
D = y(NE+4) ;

r = theta(1);
phi = theta(2);
tau = theta(3);
beta = theta(4);

etaeff = ((NE+1)/tau);

hill_factor = 1/(1+(D/Dc)^2);

dotS = r*S - phi*V*S*hill_factor ;
dotE1 = phi*S*V*hill_factor - etaeff * E_mat(1);
dotE_mat = etaeff.*E_mat(1:end-1) - etaeff.*E_mat(2:end);
dotI = etaeff * (E_mat(end) - I);
dotV = beta * etaeff * I - hill_factor*V*phi*(S+I+sum(E_mat));
dotD =  etaeff * I;

dydt = [dotS,dotE1,dotE_mat,dotI,dotV,dotD]';

end


