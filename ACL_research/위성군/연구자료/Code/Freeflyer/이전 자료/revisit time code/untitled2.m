datatmp0 = importdata(sprintf('0pointrevisitarray4.txt'))
datatmp1 = importdata(sprintf('1pointrevisitarray4.txt'))

%%
avgRT_P4_0 = mean(datatmp0.data(:,2)); 
minRT_P4_0 = min(datatmp0.data(:,2)); 
maxRT_P4_0 = max(datatmp0.data(:,2));  % Unit: day

%%
avgRT_P4_1 = mean(datatmp1.data(:,2)); 
minRT_P4_1 = min(datatmp1.data(:,2)); 
maxRT_P4_1 = max(datatmp1.data(:,2));  % Unit: day



%%
figure(); 
set(gca,'fontsize',12); hold on; grid on;
% bar([avgRT_P4*24*60, avgRT_P5*24*60, avgRT_P10*24*60]); % Unit: day-->min

bar([avgRT_P4_0*24*60, avgRT_P4_1*24*60]);
ylabel('min');
xticks([1 2 3]); 