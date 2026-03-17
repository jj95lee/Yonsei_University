function [ q_err,angle_err ] = Cal_AtdErr_2(p_i_d,aux_i_d,p_b,aux_b,q)
        
        p_b_d = quat2dcm(q)*(p_i_d);
        q_d = dcm2quat(Triad( p_i_d,aux_i_d,p_b,aux_b ));
%         q_err = qnorm(qmult([q_d(1,1) -q_d(1,2:4)],q));
		q_err = qnorm(qmult(qinv(q_d),q));
%         q_err = qnorm(qmult(q_d,q_inv(q)'));
        angle_err = acos(dot(p_b_d,p_b));
        
end

