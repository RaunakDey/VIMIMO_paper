function theta_optimized = search_minimum_onestep(theta_initial,sd_theta,data, S0,V0,num_steps, step_size)
%SEARCH_MINIMUM Summary of this function goes here
%   Detailed explanation goes here



theta = theta_initial;

for i = 1:num_steps

for j = 1:length(theta)

    theta1 = theta;
    theta2 = theta;
    
    theta1(j) = theta1(j) + sd_theta(j);
    theta2(j) = theta2(j) - sd_theta(j);
  

    error1 = error_replicates(theta1,data,S0,V0);
    error2 = error_replicates(theta2,data,S0,V0);

    

    %guaranteed minimization algorithm.

    
    [min_val,index] = min([error1,error2]);
    sign = (-1)^index;
    
    cost = error_replciates(theta,data,S0,V0)
    if min_val < cost
    theta(j) = theta(j) - sign*step_size*sd_theta(j);
    end

    %theta(j) = theta(j) + step_size*(error2 - error1)/(2*sd_theta(j));
   
end

end

theta_optimized = theta;

end

