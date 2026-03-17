function [ q_err,angle_err ] = Cal_AtdErr_1(p_i_d,p_b,q)

        p_b_d = quat2dcm(q)*(p_i_d);
        
        l = cross(p_b_d,p_b)/norm(cross(p_b_d,p_b));
        angle_err = acos(dot(p_b_d,p_b));
        q_err = [cos(angle_err/2) l(1)*sin(angle_err/2) l(2)*sin(angle_err/2) l(3)*sin(angle_err/2)];

end

