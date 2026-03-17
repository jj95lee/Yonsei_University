function qout = qinv(qin)
% q(1)   = scalar
% q(2:4) = vector
%------------------------------------------------------------------------------

[a,b] = size(qin);

if a == 4 && b ~= 4
    flag = 0;    
elseif a ~= 4 && b == 4
    flag = 1;    
else
	filename = mfilename;
    msg = ['# ',filename,' :: Dimension error'];
    error(msg);
end

qout = [qin(1);-qin(2);-qin(3);-qin(4)];
if flag
    qout = qout';
end