clear
clc
%% data setting
Timedata = importdata("f3_100sec_2days_i40_AccessReport.txt",'');
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
[a,b,c,d,e,f,g,h] = C{i,4:11};
hour = extractBetween(d,1,2);
minute = extractBetween(d,4,5);
second = extractBetween(d,7,12);
entrytime(i) = datetime(str2double(c),str2double(a),str2double(b),str2double(hour{1}),str2double(minute{1}),str2double(second{1}));
jd_entrytime(i) = juliandate(entrytime(i));
hour = extractBetween(h,1,2);
minute = extractBetween(h,4,5);
second = extractBetween(h,7,12);
exittime(i) = datetime(str2double(g),str2double(e),str2double(f),str2double(hour{1}),str2double(minute{1}),str2double(second{1})); 
jd_exittime(i) = juliandate(exittime(i));
end

%% Calculate Access Satellite Number
elapsesecond = 86400*3;
[M,I] = max(exittime);
jdtimespan = linspace(jd_entrytime(1),jd_exittime(I),elapsesecond);
K = zeros(length(jdtimespan),1);

no_sat = [];
for j = 1:length(Timedata)-1
    for i  = 1 : 1 : elapsesecond
        if jdtimespan(i) > jd_entrytime(j) &&  jdtimespan(i) < jd_exittime(j)
            K(i) = K(i)+1;
        end
    end
end

%% Calculate Discontinuity time
[sorted_entrytime,I] = sort(entrytime);
sorted_exittime = exittime(I);
sorted_time(:,1) = transpose(sorted_entrytime) ;
sorted_time(:,2) = transpose(sorted_exittime) ;
starttime = [sorted_time(1,1)];
endtime = [sorted_time(1,2)];
count = 1;

for i = 1:length(Timedata)-2
    if sorted_time(i,2) < sorted_time(i+1,1)
        count = count+1;
        starttime(count) = sorted_time(i+1,1);
        endtime(count) = sorted_time(i+1,2);
    elseif sorted_time(i,2) > sorted_time(i+1,1)
        endtime(count) = sorted_time(i+1,2);
    end
end
connect_time(:,1) = transpose(starttime) ;
connect_time(:,2) = transpose(endtime) ;
for i = 1:length(connect_time)-1
    disconnect_time(i) = connect_time(i+1,1)-connect_time(i,2);
end

A = disconnect_time
S = sum(A) %disconnect_time 총합(minute)
D = (double(minutes(S))/1440) %disconnect_time 총합(day)
P = (double(minutes(S))/4320)*100  % S값/전체 시뮬레이션 시간에 대한 % (minute)

B3 = transpose(floor(datenum(A)*86400));

Max3 = max(B3)/60
Av3 = mean(B3)/60
min3 = min(B3)/60

disconnect_time = transpose(disconnect_time);
[maxdisconnect,posi]= max(disconnect_time);
fprintf('최대 불연속 시간은 %s 입니다.\n',maxdisconnect)

disconnect_time = transpose(disconnect_time);
meandisconnect = mean(disconnect_time);
fprintf('평균 불연속 시간은 %s 입니다.\n',meandisconnect)

disconnect_time = transpose(disconnect_time);
[mindisconnect,posi]= min(disconnect_time);
fprintf('최소 불연속 시간은 %s 입니다.\n',mindisconnect)

% clearvars -except connect_time disconnect_time K jdtimespan




