clear all;
clc;



load('./CBA18-2_18.mat');
plot(time_free_phages,mean(free_phages,2),'-ko'); hold on;set(gca, 'YScale', 'log');

load('./CBA18-3_4.mat');
plot(time_free_phages,mean(free_phages,2),'-bo');

load('./CBA18-3_18.mat');
plot(time_free_phages,mean(free_phages,2),'-ro');

load('./CBA38-1_38.mat');
plot(time_free_phages,mean(free_phages,2),'-go');

load('./HP1_13-15.mat');
plot(time_free_phages,mean(free_phages,2), '-k*');

load('./HP1_H100.mat');
plot(time_free_phages,mean(free_phages,2),'-b*');

load('./HS6_13-15.mat');
plot(time_free_phages,mean(free_phages,2),'-r*');

load('./HS6_H100.mat');
plot(time_free_phages,mean(free_phages,2),'-g*');

legend('CBA18-2 -- 18','CBA18-3 -- 4','CBA18-3 -- 18','CBA38-1 -- 38','HP1 -- 13-15','HP1 -- H100','HP6 -- 13-15','HS6 -- H100'); 