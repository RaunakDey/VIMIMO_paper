function theta_optimized = search_minimum_community(theta_initial,sd_theta,data,model,pars2,mcmcpars,num_steps, step_size)
%SEARCH_MINIMUM Summary of this function goes here
%   Detailed explanation goes here

theta = theta_initial;

for i = 1:num_steps

for j = 1:length(theta)

    theta1 = theta;
    theta2 = theta;

    theta1(j) = theta1(j) + sd_theta(j);
    theta2(j) = theta2(j) - sd_theta(j);


    error1 = loglikefun(theta1,data,pars2,mcmcpars,model,0);
    error2 = loglikefun(theta2,data,pars2,mcmcpars,model,0);

    %guaranteed  (local) minimization algorithm 

    
    [min_val,index] = min([error1,error2]);
    sign = (-1)^index;
    
    cost = loglikefun(theta,data,pars2,mcmcpars,model,0)

    if min_val < cost
        if theta(j) - sign*step_size*sd_theta(j) > 0
            theta(j) = theta(j) - sign*step_size*sd_theta(j);
        end
    end

    %theta(j) = theta(j) + step_size*(error2 - error1)/(2*sd_theta(j));
   
end

end

theta_optimized = theta;

end

