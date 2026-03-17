function qerr = qerror_vec1(p1,p2)
% Quaternion error to match a couple of vector
% p1 -> p2: change attitude to match p1 to p2
% p1: present vector
% p2: target vector
% 
% q_target = qmult( q_now, q_err )

% Euler principal vector
l = cross(p1,p2);
l = l/norm(l);

% Euler principal angle
alp = acos(dot(p1,p2)/norm(p1)/norm(p2));
qerr = [cos(alp/2), l(1)*sin(alp/2), l(2)*sin(alp/2), l(3)*sin(alp/2)];

if qerr(1) < 0 
	qerr = -qerr;
end

end

