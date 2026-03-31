datatmp = importdata(sprintf('0pointrevisitarray.txt'))


%%
avgRT_P4 = mean(datatmp.data(:,2)); 
minRT_P4 = min(datatmp.data(:,2)); 
maxRT_P4 = max(datatmp.data(:,2));  % Unit: day


%%
figure(); 
set(gca,'fontsize',12); hold on; grid on;
% bar([avgRT_P4*24*60, avgRT_P5*24*60, avgRT_P10*24*60]); % Unit: day-->min

bar([avgRT_P4*24*60, minRT_P4*24*60, maxRT_P4*24*60]);
ylabel('min');
xticks([1 2 3]); 
