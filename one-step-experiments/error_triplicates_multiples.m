function ss = error_triplicates_multiples(theta,data,S0_replicates,V0_replicates)
%Error function for mcmc
% now S0 and V0 are also vectors depending on how may sets of triplicate data there is.

dilution_factor = 100;

num_replicates = length(V0_replicates)/3;


ss = 0;

for j = 1: length(V0_replicates)

S0 = S0_replicates(j);
V0 = V0_replicates(j);


% to choose which dataset
if mod(j,3) ==0
    index = 3;
else
    index = mod(j,3);
end


time_free_phages = data.xdata{ceil(j/3)}./60;

clear y0;
clear NE;
NE = round(theta(5));
y0(1) = S0;
y0(2:NE+2) = 0;
y0(NE+3) = V0;

% simulated
[~,y_series2] = one_step_simulate_particular_points(time_free_phages,y0,theta,NE,dilution_factor);
virus_simu = y_series2(end,:);

% need to put the negative numbers to 0;
virus_simu(virus_simu<0) = 0;


free_phages =  data.ydata{ceil(j/3)}(:,index);

%this is per point
% need this scaling as replicates contains unequal number of points.

%size(time_free_phages)
%size(free_phages)
%size(virus_simu')
error_this_replicate =  sum( (log(free_phages)./log(10) - log(virus_simu')./log(10) ).^2  )./length(virus_simu) ;

ss = ss + error_this_replicate;

end


ss = ss/3/num_replicates;

end



