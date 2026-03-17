function vout = qvmult(q,v)
% q(1:3) = vector
% q(4)   = scalar
%------------------------------------------------------------------------------

[a1,b1] = size(q);
[a2,b2] = size(v);

if (a1 == a2+1 && b1 == b2) || (a1 == a2 && b1 == b2+1)
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

qq = [q(1);q(2);q(3);q(4)];
vmag = norm(v);

if vmag > 0
    qv = [0;v(1);v(2);v(3)]/vmag;
    vout = qmult(qmult(qinv(qq),qv),qq);
    vout = vout(2:4,1)*vmag; 
else
    vout = zeros(3,1);
end

if flag
    vout = vout';
end