close all

figure()
set(gcf,'Color','w')
subplot(5, 1, 1)
plot(out.t_mea, out.Pgen, 'b', 'LineWidth', 1.2)
grid on
% ylim([0 100])
ylabel('Power Generation (W)')

ap = subplot(5, 1, 2)
plot(out.t_mea, out.Pload1, 'LineWidth', 1.2); hold on
plot(out.t_mea, out.Pload2, 'LineWidth', 1.2); hold on
plot(out.t_mea, out.Pload3, 'LineWidth', 1.2); hold on
plot(out.t_mea, out.Pload4, 'LineWidth', 1.2); hold on
plot(out.t_mea, out.Pload5, 'LineWidth', 1.2); hold on
plot(out.t_mea, out.Pload6, 'LineWidth', 1.2); hold on
plot(out.t_mea, out.Pload7, 'LineWidth', 1.2); hold on
plot(out.t_mea, out.Pload8, 'LineWidth', 1.2); hold on
plot(out.t_mea, out.Pload9, 'LineWidth', 1.2); hold on
plot(out.t_mea, out.Pload10, 'LineWidth', 1.2); hold on
plot(out.t_mea, out.Pload11, 'LineWidth', 1.2); hold on
plot(out.t_mea, out.Pload12, 'LineWidth', 1.2); hold on
plot(out.t_mea, out.Pload13, 'LineWidth', 1.2); hold on
plot(out.t_mea, out.Pload14, 'LineWidth', 1.2); hold on
plot(out.t_mea, out.Pload15, 'LineWidth', 1.2); hold on
plot(out.t_mea, out.Pload16, 'LineWidth', 1.2); hold on
plot(out.t_mea, out.Pload17, 'LineWidth', 1.2); hold on
plot(out.t_mea, out.Pload18, 'LineWidth', 1.2); hold on
plot(out.t_mea, out.Pload19, 'LineWidth', 1.2); hold on
plot(out.t_mea, out.Pload20, 'LineWidth', 1.2); hold 
grid on
ylim([0 7.5])
legend('EPS Battery','CSS(With Solar Panel)','CSS','CubeComputer','CubeSense','CubeControl','CubeMag','Wheel','CubeToquers(rod)','CubeToquers(coil)','GPS antenna','GPS Rx','CNDH OBC','UHF antenna','UHF Rx','UHF Tx','S-band Tx','PAYS','PAYR','PAYC')
ylabel('Power Consumption (W)')

subplot(5, 1, 3)
plot(out.t_mea, out.BatSoC, 'b', 'LineWidth', 1.2)
grid on
ylim([40 105])
ylabel('State of Charge (%)')

subplot(5, 1, 4)
plot(out.t_mea, 100 - out.BatSoC, 'b', 'LineWidth', 1.2)
grid on
ylim([0 100])
ylabel('Depth of Discharge (%)')

subplot(5, 1, 5)
plot(out.t_mea, out.Vbat, 'b', 'LineWidth', 1.2)
grid on
ylim([7 9])
ylabel('Battery Voltage (V)')
xlabel('Time (sec)')


DoD = max(out.BatSoC) - min(out.BatSoC)
minSoC = min(out.BatSoC)
maxSoC = max(out.BatSoC)
endSoC = out.BatSoC(end)

%% Power Generation
A = nnz(Pmp(1:end,1));
B = 1.135*20*A/3600*0.8*0.9*0.9; 
C = 1.135*(6*3+6*3+4*3+4*3)/(3+3+3+3+1+1)*A/3600*0.8*0.9*0.9;

if max(Pmp(:,1)) > 10
fprintf('생산전력 = %f Wh\n',B)
else 
    fprintf('생산전력 = %f Wh\n',C)
end
%% Power Consumption 
if max(Pmp(:,1)) > 10
D = (B-0.4*(out.BatSoC(end)-out.BatSoC(1)));
fprintf('소비전력 = %f Wh\n',D)
else
D = (C-0.4*(out.BatSoC(end)-out.BatSoC(1)));
fprintf('소비전력 = %f Wh\n',D)
end

%% Charging rate
E = out.BatSoC(end)-out.BatSoC(1);
fprintf('Charging rate = %f \n',E)