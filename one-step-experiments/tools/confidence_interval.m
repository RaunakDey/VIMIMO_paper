function [time,virus_all] = confidence_interval(time_free_phages,data,chain,dilution_factor)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


for i=1:length(chain)


NE = round(chain(i,5));

y0(1) = 1e8;
y0(2:NE+2) = 0;
y0(NE+3) = 1e7;

if (error_function_NE_varies(chain(i,:),data,1e8,1e7) < 2)

[time,y_series] = one_step_simulate(time_free_phages,y0,chain(i,:),NE,dilution_factor);
clear y0
clear NE
virus_all = y_series(end,:);
plot(time,virus_all,"Color",[0.7 0.7 0.7]);hold on;
set(gca,'YScale','log');
clear virus_all;


end



end

end

