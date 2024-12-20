function er_total = error_ineffi_seiv(theta,pars,model,tvec, host1, host2,host3,host4,host5,virus1,virus2,virus3,virus4,virus5)
%ERROR_INEFFI_SEIV Summary of this function goes here
%   Detailed explanation goes here



pars1 = pars;




pars1.phi = [ 0   theta(1)        0         0         0;
             theta(2)   theta(3)   theta(4)         0         0;
                0         0    theta(5)        0         0;
                0         0         0    theta(6)  theta(7);
                0         0         0   theta(8)    theta(9)];



pars1.r = [theta(10),theta(11),theta(12),theta(13),theta(14)]';



% rewrite pars1
pars1.beta =[0    theta(15)        0         0         0;
  theta(16)   theta(17)   theta(18)       0         0;
         0         0   theta(19)        0         0;
         0         0         0   theta(20)  theta(21);
         0         0         0   theta(22)  theta(23)];



pars1.tau = [0    theta(24)        0         0         0;
  theta(25)   theta(26)   theta(27)       0         0;
         0         0   theta(28)        0         0;
         0         0         0   theta(29)  theta(30);
         0         0         0   theta(31)  theta(32)];



pars1.eta(pars1.tau>0) = 1./pars1.tau(pars1.tau>0);


pars1.prob_infec  = [0    theta(33)        0         0         0;
  theta(34)   theta(35)   theta(36)       0         0;
         0         0   theta(37)        0         0;
         0         0         0   theta(38)  theta(39);
         0         0         0   theta(40)  theta(41)];


pars = pars1;


[t1,S1,V1,~] = simulate_ode_noninfec(model,pars1,tvec,pars1.S0,pars1.V0); % initial parameter set

S1(S1<1e-1) = 0;
V1(V1<1e-1) = 0;
er_host = log(S1) - log([mean(host1');mean(host2');mean(host3');mean(host4');mean(host5')]'*1e3);
er_host_t = sum(abs(er_host(~isinf(er_host))));

er_virus = log(S1) - log([mean(virus1');mean(virus2');mean(virus3');mean(virus4');mean(virus5')]'*1e3);
er_virus_t = sum(abs(er_virus(~isinf(er_virus))));

er_total = er_host_t + er_virus_t;







end

