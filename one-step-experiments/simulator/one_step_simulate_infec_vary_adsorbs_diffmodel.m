function [time,y,time_abs,pre_dil] = one_step_simulate_infec_vary_adsorbs_diffmodel(t,y0,theta,NE,dilution_factor)

opts=odeset('AbsTol',1E-8,'RelTol',1E-8,'NonNegative',1,'MaxStep',0.25);

dt_dilution = 1/60;  
time_for_absorbtion = 0:dt_dilution:15/60;  %15 mins given to us


solved = ode45(@(t,y) one_step_eqn_before_dilution(t,y,theta,NE), time_for_absorbtion',y0,opts);

%solved = ode45(@(t,y)one_step_eqn_after_dilution_infectionvary(t,y,theta,NE), time_for_absorbtion',y0,opts);

pre_dil = solved.y;
time_abs = solved.x;
y_ini = solved.y(:,end)/dilution_factor;
y_ini(end+1) = 0; % all viruses are infectitious at first.

time = 0:(t(2)-t(1))*0.01:t(end);
solved2 = ode45(@(t,y) one_step_eqn_after_dilution_infectionvary(t,y,theta,NE), time',y_ini);
y=solved2.y;
time=solved2.x;

end
