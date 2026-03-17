function qerr = qerror(q1,q2)
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

% q1 -> q2
% q1 -> ref -> q2
% q2 - q1
% q1^-1 * q2

qq1 = [q1(1);q1(2);q1(3);q1(4)];
qq2 = [q2(1);q2(2);q2(3);q2(4)];

qerr = qmult(qinv(qq1),qq2);
qerr = qnorm(qerr);

if flag
    qerr = qerr';
end