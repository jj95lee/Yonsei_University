%% Calculate Discontinuity time
startdata = importdata("start.txt",'');
enddata = importdata("end.txt",'');

count = 1;


for i = 1:2021
    if startdata(i,1) < enddata(i+1,1)
        count = count+1;
        starttime(count) = sorted_time(i+1,1);
        endtime(count) = sorted_time(i+1,2);
    elseif sorted_time(i,2) > sorted_time(i+1,1);
        endtime(count) = sorted_time(i+1,2);
    end
end
connect_time(:,1) = transpose(starttime) ;
connect_time(:,2) = transpose(endtime) ;
for i = 1:length(connect_time)-1
    disconnect_time(i) = connect_time(i+1,1)-connect_time(i,2);
end
