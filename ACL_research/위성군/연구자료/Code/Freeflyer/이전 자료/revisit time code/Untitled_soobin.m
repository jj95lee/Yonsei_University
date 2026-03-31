
% a = readtable('./revisit_1016_평면4개_1sec_섭동반영./0.0000000000000000000000000000000000000.txt');
% b = importdata('./revisit_1016_평면4개_1sec_섭동반영./0.0000000000000000000000000000000000000.txt');
% c = load('./revisit_1016_평면4개_1sec_섭동반영./0.0000000000000000000000000000000000000.txt');


%%
k = 1;
for j1 = [0 1 10 100]
    datatmp = importdata(sprintf('%d.0000000000000000000000000000000000000.txt',j1));
    dataRT_P4(k,:) = datatmp.data; k = k+1;
end


%%
for j1 = 1:86400
    datatmp = importdata(sprintf('%d.0000000000000000000000000000000000000.txt',j1));
    dataRT_P4(j1,:) = datatmp.data;
end


%%
avgRT_P4 = mean(dataRT_P4(:,2)); minRT_P4 = min(dataRT_P4(:,2)); maxRT_P4 = max(dataRT_P4(:,2));  % Unit: day


%%
figure(); 
set(gca,'fontsize',12); hold on; grid on;
% bar([avgRT_P4*24*60, avgRT_P5*24*60, avgRT_P10*24*60]);   % Unit:
% day-->min
bar([1, 2, 3]);
ylabel('min');
xticks([1 2 3]); 




