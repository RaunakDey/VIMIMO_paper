function [time,y] = one_step_simulate_seivd(t,y0,theta,NE,Dc)

opts=odeset('AbsTol',1E-8,'RelTol',1E-8,'NonNegative',1,'MaxStep',0.25);

%dt_dilution = 1/60;  
%time_for_absorbtion = 0:dt_dilution:15/60;  %15 mins given to us

% initial dead cells 
y0(end+1) = 0;

solved = ode45(@(t,y) pairwise_seivd(t,y,theta,NE,Dc), t,y0,opts);


y=solved.y;
time=solved.x;

end
