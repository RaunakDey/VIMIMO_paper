function [time,y] = simulate_pairwise_seiv(t,B0,V0,theta,NE)

opts=odeset('AbsTol',1E-8,'RelTol',1E-8,'NonNegative',1,'MaxStep',0.25);

y0(1) = B0;
y0(2:NE+2) = 0;
y0(NE+3) = V0;

length(y0)
solved = ode45(@(t,y) pairwise_seiv(t,y,theta,NE), t,y0,opts);


y=solved.y;
time=solved.x;

end
