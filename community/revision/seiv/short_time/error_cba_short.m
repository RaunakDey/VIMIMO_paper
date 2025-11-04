function error_summed = error_cba_short(theta,pars,data,model,mcmcpars,limit)
%COmputed errors from input parameters


% grab qpcr data
tvec = data.xdata;
ydata = [data.ydata(:,1:3) data.ydata(:,6:8)];

pars = update_pars(pars,theta,mcmcpars);
pars.eta(pars.tau>0) = 1./pars.tau(pars.tau>0);

S0 = pars.S0;
V0 = pars.V0;
[~,S,V] = simulate_ode(model,pars,tvec,S0,V0);
end_lim = round(length(S(:,1))*limit);

ymodel = [S V];
ymodel(ymodel<0) = 0; % get rid of negative densities

error_summed = sum(sum( (log(ymodel(1:end_lim,1:3)) - log(ydata(1:end_lim,1:3))).^2  )) + sum(sum(  (log(ymodel(:,4:6)) - log(ydata(:,4:6))).^2 )  );  



end

