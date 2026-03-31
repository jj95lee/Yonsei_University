%% Plotting
figure(1); 
set(gca,'fontsize',12); hold on; grid on;
y = [1.76 8.73 0];
b = bar([1 2 3], y);
b.FaceColor = 'flat';
b.CData(1,:) = [0 0.4470 0.7410];
b.CData(2,:) = [0.8500 0.3250 0.0980];
b.CData(3,:) = [0.9290 0.6940 0.1250];

title('T/P/F = 160:32:4 ,  Inclination = 36')
ylabel('최대 재방문주기 (min)');
xlabel('재방문 주기');
xticks([1 2 3])
xticklabels({'평균', '최대', '최소' })

