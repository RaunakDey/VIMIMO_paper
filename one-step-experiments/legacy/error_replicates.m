function [ss] = error_replicates(theta,data,S0,V0)
%Error function for mcmc

time_free_phages = data.xdata;
free_phages = data.ydata;
dilution_factor = 100;


NE = round(theta(5));
y0(1) = S0;
y0(2:NE+2) = 0;
y0(NE+3) = V0;

% simulated
[~,y_series2] = one_step_simulate_particular_points(time_free_phages,y0,theta,NE,dilution_factor);
virus_simu = y_series2(end,:);

[~,size_col] = size(free_phages);
virus_replicates = repmat(virus_simu,size_col,1)';

%ss = size(free_phages);
%ss2 =  size(virus_replicates);

ss = sum(sum( (log(free_phages)./log(10) - log(virus_replicates)./log(10)).^2) );
ss = ss/size_col; %normalize for number of replicates.

end



