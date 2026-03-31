%% Plotting
figure(1); 
set(gca,'fontsize',12); hold on; grid on;
y = [2.58 2.45 0.84 0.86 1.76];
b = bar([0 1 2 3 4], y);
b.FaceColor = 'flat';
b.CData(1,:) = [0 0.4470 0.7410];
b.CData(2,:) = [0 0.4470 0.7410];
b.CData(3,:) = [0 0.4470 0.7410];
b.CData(4,:) = [0 0.4470 0.7410];
b.CData(5,:) = [0 0.4470 0.7410];

title('inclination = 36')
ylabel('평균 재방문주기 (min)');
xlabel('');
xticks([0 1 2 3 4])
xticklabels({'f=0', 'f=1', 'f=2' 'f=3'});
