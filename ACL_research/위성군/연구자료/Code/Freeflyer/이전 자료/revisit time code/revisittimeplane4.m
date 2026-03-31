close all
clear 
clc

%% 4plane data define
datatmp4f0 = importdata('1pointrevisitarray.txt')

%% f0 (avg,min,max)
avgRT_P4_f0 = mean(datatmp4f0.data(:,2)); 
minRT_P4_f0 = min(datatmp4f0.data(:,2)); 
maxRT_P4_f0 = max(datatmp4f0.data(:,2));  % Unit: day

%% minimum 6 col
min2RT_P4_f0 = datatmp4f0.data(86400,6);


%% maximum 5 col
max2RT_P4_f0 = datatmp4f0.data(86400,5);


%%
figure('Name','Average'); 
set(gca,'fontsize',10); hold on; grid on;
% bar([avgRT_P4*24*60, avgRT_P5*24*60, avgRT_P10*24*60]); % Unit: day-->min

bar([avgRT_P4_f0*24*60]);
xlabel('f number');
ylabel('min');
ylim([3, 4])

%%
figure('Name','min1'); 
set(gca,'fontsize',10); hold on; grid on;
% bar([avgRT_P4*24*60, avgRT_P5*24*60, avgRT_P10*24*60]); % Unit: day-->min
