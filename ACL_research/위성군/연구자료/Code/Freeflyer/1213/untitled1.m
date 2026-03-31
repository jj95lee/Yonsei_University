Timedata = importdata("nSC= 30_nPlanes=3_f=2_AccessReport_kep_1day.txt",'');
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