function ss = error_function_NE_varies(theta,data,S0,V0)
%Error function for mcmc

time_free_phages = data.xdata;
free_phages_mean = data.ydata;

dilution_factor = 100;


NE = round(theta(5));
y0(1) = S0;
y0(2:NE+2) = 0;
y0(NE+3) = V0;

% simulated
[~,y_series2] = one_step_simulate_particular_points(time_free_phages,y0,theta,NE,dilution_factor);
virus_simu = y_series2(end,:);

%[size_row,size_col] = size(virus_simu);
%virus_replicates = repmat(virus_simu,size_row,size_col);

ss = sum(sum( (log(free_phages_mean) - log(virus_simu)).^2 ));


end





