function y = one_step_simulate_pred(t,theta)

S0 =1e8;
V0 = 1e7;
NE = round(theta(5));
dilution_factor = 1000;
y0(1) = S0;
y0(2:NE+2) = 0;
y0(NE+3) = V0;


opts=odeset('AbsTol',1E-8,'RelTol',1E-8,'NonNegative',1,'MaxStep',0.25);

dt_dilution = 1/60;  
time_for_absorbtion = 0:dt_dilution:15/60;  %15 mins given to us


solved = ode45(@(t,y) one_step_eqn_before_dilution(t,y,theta,NE), time_for_absorbtion',y0,opts);

pre_dil = solved.y;
time_abs = solved.x;
y_ini = solved.y(:,end)/dilution_factor;
%y_ini(end) = 2e4; %forcibly putting to initial value approx
%y_ini(3:NE+1) =0; %forcibly putting to 0

time = 0:(t(2)-t(1))*0.01:t(end);
solved2 = ode45(@(t,y) one_step_eqn_before_dilution(t,y,theta,NE), time',y_ini);
y=solved2.y;
time=solved2.x;

end
