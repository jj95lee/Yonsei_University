%% Plotting
figure(1); 
set(gca,'fontsize',12); hold on; grid on;
y = [2.58 2.43 0.05 0 0];
b = bar([0 1 2 3 4], y);
b.FaceColor = 'flat';
b.CData(1,:) = [0.9290 0.6940 0.1250];
b.CData(2,:) = [0.9290 0.6940 0.1250];
b.CData(3,:) = [0.9290 0.6940 0.1250];
b.CData(4,:) = [0.9290 0.6940 0.1250];
b.CData(5,:) = [0.9290 0.6940 0.1250];

title('inclination = 36')
ylabel('최소 재방문주기 (min)');
xlabel('');
xticks([0 1 2 3 4])
xticklabels({'f=0', 'f=1', 'f=2' 'f=3'});