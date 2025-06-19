function [time,y] = one_step_simulate_noninfectious_seivd(t,B0,V0,theta,NE,Dc,infectious_ratio)

opts=odeset('AbsTol',1E-8,'RelTol',1E-8,'NonNegative',1,'MaxStep',0.25);

y0(1) = B0;
y0(2:NE+2) = 0;
y0(NE+3) = V0*infectious_ratio;
y0(NE+4) = V0*(1-infectious_ratio);
y0(NE+5) = 0; % initial dead cells 

solved = ode45(@(t,y) pairwise_noninfectious_seivd(t,y,theta,NE,Dc), t,y0,opts);


y=solved.y;
time=solved.x;

end
