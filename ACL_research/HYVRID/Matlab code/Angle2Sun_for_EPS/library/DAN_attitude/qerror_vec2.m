function qerr = qerror_vec2(p11,p12,p21,p22)
% Quaternion error to match a couple of vector
% p11 -> p21: change attitude to match p11 to p21
% p11:	present primary vector
% p12:	present secondary vector
% p21:	target primary vector
% p22:	target secondary vector
%
% q_target = qmult( q_now, q_err )

qerr = dcm2quat(triad( p11,p12,p21,p22 ));

if qerr(1) < 0 
	qerr = -qerr;
end

end

