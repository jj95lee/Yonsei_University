clear
clc
% data setting
Timedata = importdata("0_100sec_7days_AccessReport.txt",'');

k = Timedata(2:length(Timedata),1);
C = {};
entrytime = datetime();
exittime = datetime();
jd_entrytime = zeros(length(Timedata)-1,1);
jd_exittime = zeros(length(Timedata)-1,1);
for i = 1: length(Timedata)-1
    C(i,:) = strsplit(k{i});
    for j = 1:12
        switch C{i,j}
            case 'Jan'
                C{i,j} = '01';
        end
    end