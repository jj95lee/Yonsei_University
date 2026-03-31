startdata = importdata("start.txt",'');
enddata = importdata("end.txt",'');

connect_time(:,1) = transpose(startdata(1,:)) ;
connect_time(:,2) = transpose(enddata(1,:)) ;
for i = 1:length(connect_time)-1
    disconnect_time(i) = connect_time(i+1,1)-connect_time(i,2);
end