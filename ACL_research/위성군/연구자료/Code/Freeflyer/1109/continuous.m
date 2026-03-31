clear
clc
%% data setting
Timedata = importdata("nSC= 44_nPlanes=4_f=0_AccessReportnew.txt",'');
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
min = extractBetween(d,4,5);
second = extractBetween(d,7,12);
entrytime(i) = datetime(str2double(c),str2double(a),str2double(b),str2double(hour{1}),str2double(min{1}),str2double(second{1}));
jd_entrytime(i) = juliandate(entrytime(i));
hour = extractBetween(h,1,2);
min = extractBetween(h,4,5);
second = extractBetween(h,7,12);
exittime(i) = datetime(str2double(g),str2double(e),str2double(f),str2double(hour{1}),str2double(min{1}),str2double(second{1})); 
jd_exittime(i) = juliandate(exittime(i));
end

% Calculate Access Satellite Number
elapsesecond = 86400*7;
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

% Calculate Discontinuity time
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

disconnect_time = transpose(disconnect_time);
[maxdisconnect,posi]= max(disconnect_time);
fprintf('최대 불연속 시간은 %s 입니다.\n',maxdisconnect)
clearvars -except connect_time disconnect_time K jdtimespan

%Plotting
figure(1)
histogram(disconnect_time,200)

figure(2)
timeaa = datetime(jdtimespan,'ConvertFrom','juliandate');
plot(timeaa,K)
