function [ u_atd ] = Ctrl_Atd_PD( q_err,w_err,T_max,ctrlgain_atd_PD )
%#codegen
% PD control for attitude tracking with reaction wheels

kp = ctrlgain_atd_PD(1:3);
kd = ctrlgain_atd_PD(4:6);

u_atd = -kp.*q_err(2:4)' - kd.*w_err;

for j=1:3
    if abs(u_atd(j)) > T_max
        u_atd(j) = sign(u_atd(j))*T_max;
    end
end

end

