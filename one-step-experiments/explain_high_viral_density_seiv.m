clc;
clear all;
addpath(genpath(pwd)); % add current directory to path (mcmcstat is included)
addpath('./simulator/');
addpath('./mcmcstat/');


%% 1

cutoff = 1e0;
tvec = 0:0.1:15.8;
MOI_list = [10,1,0.1];
color_list = ['r','b','k'];

for iter = 1:length(MOI_list)

MOI  = MOI_list(iter);
color_line = color_list(iter);

B0 = 2e6;
V0 = B0*MOI;
load('./data_2024/CBA18-2_18_2024.mat','theta_optimized');
NE_optimal = round(theta_optimized(5));



[time_seiv, y_series_seiv] = simulate_pairwise_seiv(tvec,B0,V0,theta_optimized,NE_optimal);

virus_den_seiv = y_series_seiv(NE_optimal+3,:);
host_den_seiv = sum(y_series_seiv(1:NE_optimal+2,:));
host_den_seiv(host_den_seiv<cutoff) = 0;
  

subplot(4,4,1)
plot(time_seiv,host_den_seiv, Color=color_line)
set(gca,'YScale','log');
xlabel('Time (hr)');
ylabel('Host/ml');
set(gca,'FontSize',20);
yticks([1e-2 1e0 1e2 1e4 1e6 1e8 1e10]);
hold on;
title('CBA18 -- CBA18:2');
 
subplot(4,4,2)
plot(time_seiv,virus_den_seiv, Color=color_line);
set(gca,'YScale','log');
xlabel('Time (hr)');
ylabel('Phage/ml');
set(gca,'FontSize',20);
yticks([1e-2 1e0 1e2 1e4 1e6 1e8 1e10]);
title('CBA18 -- CBA18:2');
hold on;

end

clear all;

%% 2

cutoff = 1e-2;
tvec = 0:0.1:15.8;
MOI_list = [10,1,0.1];
color_list = ['r','b','k'];

for iter = 1:length(MOI_list)

MOI  = MOI_list(iter);
color_line = color_list(iter);

B0 = 2e6;
V0 = B0*MOI;
load('./data_2024/CBA18-3_4_2024.mat','theta_optimized');
NE_optimal = round(theta_optimized(5));



[time_seiv, y_series_seiv] = simulate_pairwise_seiv(tvec,B0,V0,theta_optimized,NE_optimal);

virus_den_seiv = y_series_seiv(NE_optimal+3,:);
host_den_seiv = sum(y_series_seiv(1:NE_optimal+2,:));
host_den_seiv(host_den_seiv<cutoff) = 0;
  

subplot(4,4,3)
plot(time_seiv,host_den_seiv, Color=color_line)
set(gca,'YScale','log');
xlabel('Time (hr)');
ylabel('Host/ml');
set(gca,'FontSize',20);
yticks([1e-2 1e0 1e2 1e4 1e6 1e8 1e10]);
title('CBA4 -- CBA18:3');
hold on;

 
subplot(4,4,4)
plot(time_seiv,virus_den_seiv, Color=color_line);
set(gca,'YScale','log');
xlabel('Time (hr)');
ylabel('Phage/ml');
set(gca,'FontSize',20);
yticks([1e-2 1e0 1e2 1e4 1e6 1e8 1e10]);
title('CBA4 -- CBA18:3');
hold on;

end

clear all;

%% 3

cutoff = 1e-2;
tvec = 0:0.1:15.8;
MOI_list = [10,1,0.1];
color_list = ['r','b','k'];

for iter = 1:length(MOI_list)

MOI  = MOI_list(iter);
color_line = color_list(iter);

B0 = 2e6;
V0 = B0*MOI;
load('./data_2024/CBA18-3_18_2024.mat','theta_optimized');
NE_optimal = round(theta_optimized(5));



[time_seiv, y_series_seiv] = simulate_pairwise_seiv(tvec,B0,V0,theta_optimized,NE_optimal);

virus_den_seiv = y_series_seiv(NE_optimal+3,:);
host_den_seiv = sum(y_series_seiv(1:NE_optimal+2,:));
host_den_seiv(host_den_seiv<cutoff) = 0;
  

subplot(4,4,5)
plot(time_seiv,host_den_seiv, Color=color_line)
set(gca,'YScale','log');
xlabel('Time (hr)');
ylabel('Host/ml');
set(gca,'FontSize',20);
yticks([1e-2 1e0 1e2 1e4 1e6 1e8 1e10]);
title('CBA18 -- CBA18:3');
hold on;

 
subplot(4,4,6)
plot(time_seiv,virus_den_seiv, Color=color_line);
set(gca,'YScale','log');
xlabel('Time (hr)');
ylabel('Phage/ml');
set(gca,'FontSize',20);
yticks([1e-2 1e0 1e2 1e4 1e6 1e8 1e10]);
title('CBA18 -- CBA18:3');
hold on;

end

clear all;

%% 4

cutoff = 1e-2;
tvec = 0:0.1:15.8;
MOI_list = [10,1,0.1];
color_list = ['r','b','k'];

for iter = 1:length(MOI_list)

MOI  = MOI_list(iter);
color_line = color_list(iter);

B0 = 2e6;
V0 = B0*MOI;
load('./data_2024/CBA38-1_38.mat','theta_optimized');
NE_optimal = round(theta_optimized(5));



[time_seiv, y_series_seiv] = simulate_pairwise_seiv(tvec,B0,V0,theta_optimized,NE_optimal);

virus_den_seiv = y_series_seiv(NE_optimal+3,:);
host_den_seiv = sum(y_series_seiv(1:NE_optimal+2,:));
host_den_seiv(host_den_seiv<cutoff) = 0;
  

subplot(4,4,7)
plot(time_seiv,host_den_seiv, Color=color_line)
set(gca,'YScale','log');
xlabel('Time (hr)');
ylabel('Host/ml');
set(gca,'FontSize',20);
yticks([1e-2 1e0 1e2 1e4 1e6 1e8 1e10]);
title('CBA38 -- CBA38:1');
hold on;

 
subplot(4,4,8)
plot(time_seiv,virus_den_seiv, Color=color_line);
set(gca,'YScale','log');
xlabel('Time (hr)');
ylabel('Phage/ml');
set(gca,'FontSize',20);
yticks([1e-2 1e0 1e2 1e4 1e6 1e8 1e10]);
title('CBA38 -- CBA38:1');
hold on;

end

clear all;

%% 5

cutoff = 1e-2;
tvec = 0:0.1:15.8;
MOI_list = [10,1,0.1];
color_list = ['r','b','k'];

for iter = 1:length(MOI_list)

MOI  = MOI_list(iter);
color_line = color_list(iter);

B0 = 2e6;
V0 = B0*MOI;
load('./data_2024/HP1_H100_2024.mat','theta_optimized');
NE_optimal = round(theta_optimized(5));



[time_seiv, y_series_seiv] = simulate_pairwise_seiv(tvec,B0,V0,theta_optimized,NE_optimal);

virus_den_seiv = y_series_seiv(NE_optimal+3,:);
host_den_seiv = sum(y_series_seiv(1:NE_optimal+2,:));
host_den_seiv(host_den_seiv<cutoff) = 0;
  

subplot(4,4,9)
plot(time_seiv,host_den_seiv, Color=color_line)
set(gca,'YScale','log');
xlabel('Time (hr)');
ylabel('Host/ml');
set(gca,'FontSize',20);
yticks([1e-2 1e0 1e2 1e4 1e6 1e8 1e10]);
title('HP1 -- H100');
hold on;

 
subplot(4,4,10)
plot(time_seiv,virus_den_seiv, Color=color_line);
set(gca,'YScale','log');
xlabel('Time (hr)');
ylabel('Phage/ml');
set(gca,'FontSize',20);
yticks([1e-2 1e0 1e2 1e4 1e6 1e8 1e10]);
title('HP1 -- H100');
hold on;

end

clear all;

%% 6



cutoff = 1e-2;
tvec = 0:0.1:15.8;
MOI_list = [10,1,0.1];
color_list = ['r','b','k'];

for iter = 1:length(MOI_list)

MOI  = MOI_list(iter);
color_line = color_list(iter);

B0 = 2e6;
V0 = B0*MOI;
load('./data_2024/HP1_13-15_2024.mat','theta_optimized');
NE_optimal = round(theta_optimized(5));



[time_seiv, y_series_seiv] = simulate_pairwise_seiv(tvec,B0,V0,theta_optimized,NE_optimal);

virus_den_seiv = y_series_seiv(NE_optimal+3,:);
host_den_seiv = sum(y_series_seiv(1:NE_optimal+2,:));
host_den_seiv(host_den_seiv<cutoff) = 0;
  

subplot(4,4,11)
plot(time_seiv,host_den_seiv, Color=color_line)
set(gca,'YScale','log');
xlabel('Time (hr)');
ylabel('Host/ml');
set(gca,'FontSize',20);
yticks([1e-2 1e0 1e2 1e4 1e6 1e8 1e10]);
title('HP1 -- 13-15');
hold on;

 
subplot(4,4,12)
plot(time_seiv,virus_den_seiv, Color=color_line);
set(gca,'YScale','log');
xlabel('Time (hr)');
ylabel('Phage/ml');
set(gca,'FontSize',20);
yticks([1e-2 1e0 1e2 1e4 1e6 1e8 1e10]);
title('HP1 -- 13-15');
hold on;

end

clear all;

%% 7



cutoff = 1e-2;
tvec = 0:0.1:15.8;
MOI_list = [10,1,0.1];
color_list = ['r','b','k'];

for iter = 1:length(MOI_list)

MOI  = MOI_list(iter);
color_line = color_list(iter);

B0 = 2e6;
V0 = B0*MOI;
load('./data_2024/HS6_H100_2024.mat','theta_optimized');
NE_optimal = round(theta_optimized(5));



[time_seiv, y_series_seiv] = simulate_pairwise_seiv(tvec,B0,V0,theta_optimized,NE_optimal);

virus_den_seiv = y_series_seiv(NE_optimal+3,:);
host_den_seiv = sum(y_series_seiv(1:NE_optimal+2,:));
host_den_seiv(host_den_seiv<cutoff) = 0;
  

subplot(4,4,13)
plot(time_seiv,host_den_seiv, Color=color_line)
set(gca,'YScale','log');
xlabel('Time (hr)');
ylabel('Host/ml');
set(gca,'FontSize',20);
yticks([1e-2 1e0 1e2 1e4 1e6 1e8 1e10]);
title('HS6 -- H100');
hold on;

 
subplot(4,4,14)
plot(time_seiv,virus_den_seiv, Color=color_line);
set(gca,'YScale','log');
xlabel('Time (hr)');
ylabel('Phage/ml');
set(gca,'FontSize',20);
yticks([1e-2 1e0 1e2 1e4 1e6 1e8 1e10]);
title('HS6 -- H100');
hold on;

end

clear all;


%% 8

cutoff = 1e-2;
tvec = 0:0.1:15.8;
MOI_list = [10,1,0.1];
color_list = ['r','b','k'];

for iter = 1:length(MOI_list)

MOI  = MOI_list(iter);
color_line = color_list(iter);

B0 = 2e6;
V0 = B0*MOI;
load('./data_2024/HS6_13-15_2024.mat','theta_optimized');
NE_optimal = round(theta_optimized(5));



[time_seiv, y_series_seiv] = simulate_pairwise_seiv(tvec,B0,V0,theta_optimized,NE_optimal);

virus_den_seiv = y_series_seiv(NE_optimal+3,:);
host_den_seiv = sum(y_series_seiv(1:NE_optimal+2,:));
host_den_seiv(host_den_seiv<cutoff) = 0;
  

subplot(4,4,15)
plot(time_seiv,host_den_seiv, Color=color_line)
set(gca,'YScale','log');
xlabel('Time (hr)');
ylabel('Host/ml');
set(gca,'FontSize',20);
yticks([1e-2 1e0 1e2 1e4 1e6 1e8 1e10]);
title('HS6 -- 13-15');
hold on;

 
subplot(4,4,16)
plot(time_seiv,virus_den_seiv, Color=color_line);
set(gca,'YScale','log');
xlabel('Time (hr)');
ylabel('Phage/ml');
set(gca,'FontSize',20);
yticks([1e-2 1e0 1e2 1e4 1e6 1e8 1e10]);
title('HS6 -- 13-15');
hold on;

end

clear all;



%% touch up

%% 1 

cutoff = 1e-2;
tvec = 0:0.1:15.8;
MOI_list = [10,1,0.1];
color_list = ['r','b','k'];

for iter = 1:length(MOI_list)

MOI  = MOI_list(iter);
color_line = color_list(iter);

B0 = 2e6;
V0 = B0*MOI;
load('./data_2024/CBA18-2_18_2024.mat','theta_optimized');
NE_optimal = round(theta_optimized(5));
Dc = 6.19e6;
Dc = 1e5;


[time_seivd, y_series_seivd] = simulate_pairwise_seivd(tvec,B0,V0,theta_optimized,NE_optimal,Dc);

virus_den_seivd = y_series_seivd(NE_optimal+3,:);
host_den_seivd = sum(y_series_seivd(1:NE_optimal+2,:));
host_den_seivd(host_den_seivd<cutoff) = 0;
  

subplot(4,4,1)
plot(time_seivd,host_den_seivd, Color=color_line,LineStyle='--');
set(gca,'YScale','log');
xlabel('Time (hr)');
ylabel('Host/ml');
set(gca,'FontSize',20);
yticks([1e-2 1e0 1e2 1e4 1e6 1e8 1e10]);
hold on;

 
subplot(4,4,2)
plot(time_seivd,virus_den_seivd, Color=color_line,LineStyle='--')
set(gca,'YScale','log');
xlabel('Time (hr)');
ylabel('Phage/ml');
set(gca,'FontSize',20);
yticks([1e-2 1e0 1e2 1e4 1e6 1e8 1e10]);
hold on;

end

clear all;


%% 2

cutoff = 1e-2;
tvec = 0:0.1:15.8;
MOI_list = [10,1,0.1];
color_list = ['r','b','k'];

for iter = 1:length(MOI_list)

MOI  = MOI_list(iter);
color_line = color_list(iter);

B0 = 2e6;
V0 = B0*MOI;
load('./data_2024/CBA18-3_4_2024.mat','theta_optimized');
NE_optimal = round(theta_optimized(5));
Dc = 5e6;
Dc = 1e5;


[time_seivd, y_series_seivd] = simulate_pairwise_seivd(tvec,B0,V0,theta_optimized,NE_optimal,Dc);

virus_den_seivd = y_series_seivd(NE_optimal+3,:);
host_den_seivd = sum(y_series_seivd(1:NE_optimal+2,:));
host_den_seivd(host_den_seivd<cutoff) = 0;
  

subplot(4,4,3)
plot(time_seivd,host_den_seivd, Color=color_line,LineStyle='--')
set(gca,'YScale','log');
xlabel('Time (hr)');
ylabel('Host/ml');
set(gca,'FontSize',20);
yticks([1e-2 1e0 1e2 1e4 1e6 1e8 1e10]);
hold on;

 
subplot(4,4,4)
plot(time_seivd,virus_den_seivd, Color=color_line,LineStyle='--')
set(gca,'YScale','log');
xlabel('Time (hr)');
ylabel('Phage/ml');
set(gca,'FontSize',20);
yticks([1e-2 1e0 1e2 1e4 1e6 1e8 1e10]);
hold on;

end

clear all;


%% 3


cutoff = 1e-2;
tvec = 0:0.1:15.8;
MOI_list = [10,1,0.1];
color_list = ['r','b','k'];

for iter = 1:length(MOI_list)

MOI  = MOI_list(iter);
color_line = color_list(iter);

B0 = 2e6;
V0 = B0*MOI;
load('./data_2024/CBA18-3_18_2024.mat','theta_optimized');
NE_optimal = round(theta_optimized(5));
Dc = 6.19e6;
Dc = 1e5;


[time_seivd, y_series_seivd] = simulate_pairwise_seivd(tvec,B0,V0,theta_optimized,NE_optimal,Dc);

virus_den_seivd = y_series_seivd(NE_optimal+3,:);
host_den_seivd = sum(y_series_seivd(1:NE_optimal+2,:));
host_den_seivd(host_den_seivd<cutoff) = 0;
  

subplot(4,4,5)
plot(time_seivd,host_den_seivd, Color=color_line,LineStyle='--')
set(gca,'YScale','log');
xlabel('Time (hr)');
ylabel('Host/ml');
set(gca,'FontSize',20);
yticks([1e-2 1e0 1e2 1e4 1e6 1e8 1e10]);
hold on;

 
subplot(4,4,6)
plot(time_seivd,virus_den_seivd, Color=color_line,LineStyle='--')
set(gca,'YScale','log');
xlabel('Time (hr)');
ylabel('Phage/ml');
set(gca,'FontSize',20);
yticks([1e-2 1e0 1e2 1e4 1e6 1e8 1e10]);
hold on;

end

clear all;



%% 4


cutoff = 1e-2;
tvec = 0:0.1:15.8;
MOI_list = [10,1,0.1];
color_list = ['r','b','k'];

for iter = 1:length(MOI_list)

MOI  = MOI_list(iter);
color_line = color_list(iter);

B0 = 2e6;
V0 = B0*MOI;
load('./data_2024/CBA38-1_38.mat','theta_optimized');
NE_optimal = round(theta_optimized(5));
Dc = 11.26e6;
Dc = 1e5;


[time_seivd, y_series_seivd] = simulate_pairwise_seivd(tvec,B0,V0,theta_optimized,NE_optimal,Dc);

virus_den_seivd = y_series_seivd(NE_optimal+3,:);
host_den_seivd = sum(y_series_seivd(1:NE_optimal+2,:));
host_den_seivd(host_den_seivd<cutoff) = 0;
  

subplot(4,4,7)
plot(time_seivd,host_den_seivd, Color=color_line,LineStyle='--')
set(gca,'YScale','log');
xlabel('Time (hr)');
ylabel('Host/ml');
set(gca,'FontSize',20);
yticks([1e-2 1e0 1e2 1e4 1e6 1e8 1e10]);
hold on;

 
subplot(4,4,8)
plot(time_seivd,virus_den_seivd, Color=color_line,LineStyle='--')
set(gca,'YScale','log');
xlabel('Time (hr)');
ylabel('Phage/ml');
set(gca,'FontSize',20);
yticks([1e-2 1e0 1e2 1e4 1e6 1e8 1e10]);
hold on;

end

clear all;

%% 5


cutoff = 1e-2;
tvec = 0:0.1:15.8;
MOI_list = [10,1,0.1];
color_list = ['r','b','k'];

for iter = 1:length(MOI_list)

MOI  = MOI_list(iter);
color_line = color_list(iter);

B0 = 2e6;
V0 = B0*MOI;
load('./data_2024/HP1_H100_2024.mat','theta_optimized');
NE_optimal = round(theta_optimized(5));
Dc = 2.16e6;
Dc = 1e5;


[time_seivd, y_series_seivd] = simulate_pairwise_seivd(tvec,B0,V0,theta_optimized,NE_optimal,Dc);

virus_den_seivd = y_series_seivd(NE_optimal+3,:);
host_den_seivd = sum(y_series_seivd(1:NE_optimal+2,:));
host_den_seivd(host_den_seivd<cutoff) = 0;
  

subplot(4,4,9)
plot(time_seivd,host_den_seivd, Color=color_line,LineStyle='--')
set(gca,'YScale','log');
xlabel('Time (hr)');
ylabel('Host/ml');
set(gca,'FontSize',20);
yticks([1e-2 1e0 1e2 1e4 1e6 1e8 1e10]);
hold on;

 
subplot(4,4,10)
plot(time_seivd,virus_den_seivd, Color=color_line,LineStyle='--')
set(gca,'YScale','log');
xlabel('Time (hr)');
ylabel('Phage/ml');
set(gca,'FontSize',20);
yticks([1e-2 1e0 1e2 1e4 1e6 1e8 1e10]);
hold on;

end

clear all;

%% 6


cutoff = 1e-2;
tvec = 0:0.1:15.8;
MOI_list = [10,1,0.1];
color_list = ['r','b','k'];

for iter = 1:length(MOI_list)

MOI  = MOI_list(iter);
color_line = color_list(iter);

B0 = 2e6;
V0 = B0*MOI;
load('./data_2024/HP1_13-15_2024.mat','theta_optimized');
NE_optimal = round(theta_optimized(5));
Dc = 1.61e6;
Dc = 1e5;


[time_seivd, y_series_seivd] = simulate_pairwise_seivd(tvec,B0,V0,theta_optimized,NE_optimal,Dc);

virus_den_seivd = y_series_seivd(NE_optimal+3,:);
host_den_seivd = sum(y_series_seivd(1:NE_optimal+2,:));
host_den_seivd(host_den_seivd<cutoff) = 0;
  

subplot(4,4,11)
plot(time_seivd,host_den_seivd, Color=color_line,LineStyle='--')
set(gca,'YScale','log');
xlabel('Time (hr)');
ylabel('Host/ml');
set(gca,'FontSize',20);
yticks([1e-2 1e0 1e2 1e4 1e6 1e8 1e10]);
hold on;

 
subplot(4,4,12)
plot(time_seivd,virus_den_seivd, Color=color_line,LineStyle='--')
set(gca,'YScale','log');
xlabel('Time (hr)');
ylabel('Phage/ml');
set(gca,'FontSize',20);
yticks([1e-2 1e0 1e2 1e4 1e6 1e8 1e10]);
hold on;

end

clear all;


%% 7 

cutoff = 1e-2;
tvec = 0:0.1:15.8;
MOI_list = [10,1,0.1];
color_list = ['r','b','k'];

for iter = 1:length(MOI_list)

MOI  = MOI_list(iter);
color_line = color_list(iter);

B0 = 2e6;
V0 = B0*MOI;
load('./data_2024/HS6_13-15_2024.mat','theta_optimized');
NE_optimal = round(theta_optimized(5));
Dc = 2.16e6;
Dc = 1e5;

[time_seivd, y_series_seivd] = simulate_pairwise_seivd(tvec,B0,V0,theta_optimized,NE_optimal,Dc);

virus_den_seivd = y_series_seivd(NE_optimal+3,:);
host_den_seivd = sum(y_series_seivd(1:NE_optimal+2,:));
host_den_seivd(host_den_seivd<cutoff) = 0;
  

subplot(4,4,13)
plot(time_seivd,host_den_seivd, Color=color_line,LineStyle='--')
set(gca,'YScale','log');
xlabel('Time (hr)');
ylabel('Host/ml');
set(gca,'FontSize',20);
yticks([1e-2 1e0 1e2 1e4 1e6 1e8 1e10]);
hold on;

 
subplot(4,4,14)
plot(time_seivd,virus_den_seivd, Color=color_line,LineStyle='--')
set(gca,'YScale','log');
xlabel('Time (hr)');
ylabel('Phage/ml');
set(gca,'FontSize',20);
yticks([1e-2 1e0 1e2 1e4 1e6 1e8 1e10]);
hold on;

end

clear all;

%% 8
cutoff = 1e-2;
tvec = 0:0.1:15.8;
MOI_list = [10,1,0.1];
color_list = ['r','b','k'];

for iter = 1:length(MOI_list)

MOI  = MOI_list(iter);
color_line = color_list(iter);

B0 = 2e6;
V0 = B0*MOI;
load('./data_2024/HS6_13-15_2024.mat','theta_optimized');
NE_optimal = round(theta_optimized(5));
Dc = 1.61e6;
Dc = 1e5;

[time_seivd, y_series_seivd] = simulate_pairwise_seivd(tvec,B0,V0,theta_optimized,NE_optimal,Dc);

virus_den_seivd = y_series_seivd(NE_optimal+3,:);
host_den_seivd = sum(y_series_seivd(1:NE_optimal+2,:));
host_den_seivd(host_den_seivd<cutoff) = 0;
  

subplot(4,4,15)
plot(time_seivd,host_den_seivd, Color=color_line,LineStyle='--')
set(gca,'YScale','log');
xlabel('Time (hr)');
ylabel('Host/ml');
set(gca,'FontSize',20);
yticks([1e-2 1e0 1e2 1e4 1e6 1e8 1e10]);
hold on;

 
subplot(4,4,16)
plot(time_seivd,virus_den_seivd, Color=color_line,LineStyle='--')
set(gca,'YScale','log');
xlabel('Time (hr)');
ylabel('Phage/ml');
set(gca,'FontSize',20);
yticks([1e-2 1e0 1e2 1e4 1e6 1e8 1e10]);
hold on;

end

clear all;

legend('seiv moi=10','seiv moi=1','seiv moi=0.1','seivd moi=10','seivd moi=1','seivd moi=0.1')