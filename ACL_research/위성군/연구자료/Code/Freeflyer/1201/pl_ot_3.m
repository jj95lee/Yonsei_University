%% Plotting
figure(1); 
set(gca,'fontsize',12); hold on; grid on;
y = [2.58 2.47 1.6 7.97 8.73];
b = bar([0 1 2 3 4], y);
b.FaceColor = 'flat';
b.CData(1,:) = [0.8500 0.3250 0.0980];
b.CData(2,:) = [0.8500 0.3250 0.0980];
b.CData(3,:) = [0.8500 0.3250 0.0980];
b.CData(4,:) = [0.8500 0.3250 0.0980];
b.CData(5,:) = [0.8500 0.3250 0.0980];

title('inclination = 36')
ylabel('최대 재방문주기 (min)');
xlabel('');
xticks([0 1 2 3 4])
xticklabels({'f=0', 'f=1', 'f=2' 'f=3'});