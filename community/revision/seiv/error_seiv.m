function error_summed = error_seiv(theta,pars,data,model,mcmcpars,limit)
%COmputed errors from input parameters


% grab qpcr data
tvec = data.xdata;
ydata = data.ydata;

pars = update_pars(pars,theta,mcmcpars);

S0 = pars.S0;
V0 = pars.V0;
[~,S,V] = simulate_ode(model,pars,tvec,S0,V0);

ymodel = [S V];
ymodel(ymodel<0) = 0; % get rid of negative densities

end_lim = round(length(S(:,1))*limit);

error_summed = sum(sum((log(ymodel(1:end_lim,1:5)) - log(ydata(1:end_lim,1:5))).^2)) + sum(sum((log(ymodel(1:end,6:10)) - log(ydata(1:end,6:10))).^2));  


end

