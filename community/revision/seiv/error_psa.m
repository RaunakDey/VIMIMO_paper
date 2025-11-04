function error_summed = error_psa(theta,pars,data,model,mcmcpars)
%COmputed errors from input parameters


% grab qpcr data
tvec = data.xdata;
ydata = [data.ydata(:,4:5) data.ydata(:,9:10)];

pars = update_pars(pars,theta,mcmcpars);
pars.eta(pars.tau>0) = 1./pars.tau(pars.tau>0);

S0 = pars.S0;
V0 = pars.V0;
[~,S,V] = simulate_ode(model,pars,tvec,S0,V0);

ymodel = [S V];
ymodel(ymodel<0) = 0; % get rid of negative densities

error_summed = sum(sum((log(ydata) - log(ymodel)).^2));



end

