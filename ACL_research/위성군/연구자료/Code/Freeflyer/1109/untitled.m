%% Call continuous
continuous_p4_f1
continuous_p4_f2
continuous_p4_f3

%% Plotting
figure(1); 
set(gca,'fontsize',12); hold on; grid on;
y = [Max0 Max1 Max2 Max3];
b = bar([0 1 2 3], y);
b.FaceColor = 'flat';
b.CData(1,:) = [0 0.4470 0.7410];
b.CData(2,:) = [0.8500 0.3250 0.0980];
b.CData(3,:) = [0.9290 0.6940 0.1250];
b.CData(4,:) = [0.4940 0.1840 0.5560];

ylabel('최대 재방문주기 (min)');
xlabel('T:P:F');
xticks([0 1 2 3])
xticklabels({'160:4:0', '160:4:1', '160:4:2', '160:4:3'})


