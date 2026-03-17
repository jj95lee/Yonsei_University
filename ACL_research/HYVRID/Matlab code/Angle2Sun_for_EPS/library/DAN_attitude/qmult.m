function qout = qmult(q1,q2)
% q(1)   = scalar
% q(2:4) = vector
%------------------------------------------------------------------------------

[a1,b1] = size(q1);
[a2,b2] = size(q2);

if a1 == a2 && b1 == b2
    if a1 == 4 && b1 ~= 4
        flag = 0;
    elseif a1 ~= 4 && b1 == 4
        flag = 1;
	else
		filename = mfilename;
        msg = ['# ',filename,' :: Dimension error'];
        error(msg);
    end
else
    filename = mfilename;
	msg = ['# ',filename,' :: Dimension error'];
    error(msg);
end

q1s = q1(1);
q1v = [q1(2);q1(3);q1(4)];
q2s = q2(1);
q2v = [q2(2);q2(3);q2(4)];

qs = q1s*q2s - q1v'*q2v;
qv = q1s*q2v + q1v*q2s + cross(q1v,q2v);
qout = [qs;qv];
qout = qnorm(qout);

if flag
    qout = qout';
end