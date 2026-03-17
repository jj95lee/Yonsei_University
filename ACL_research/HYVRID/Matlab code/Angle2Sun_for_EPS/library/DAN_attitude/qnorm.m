function qout = qnorm(qin)
% quaternion normalization code
% q(1:3) = vector
% q(4)   = scalar
%------------------------------------------------------------------------------

[a,b] = size(qin);

if a == 4 && b ~= 4
    flag = 0;    
elseif a ~= 4 && b == 4
    flag = 1;    
else
    msg = '# qnorm.m :: Dimension error';
    error(msg);
end

qout = [qin(1);qin(2);qin(3);qin(4)]/sqrt(qin(1)^2+qin(2)^2+qin(3)^2+qin(4)^2);
if qout(1) < 0
	qout = -qout;
end

if flag
    qout = qout';
end