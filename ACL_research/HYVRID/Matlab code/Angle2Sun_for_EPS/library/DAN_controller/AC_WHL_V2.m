function [T, q_err] = AC_WHL_V2( t,eX_atd,eX_orb,ref_atd,quat_targ,auxdata,scdata,flag_AC,flag_fss,submode )

MOI = scdata.MOI;
T_max = scdata.T_max;
ctrlgain_atd_PD = scdata.ctrlgain_atd_PD;

rSun_i = ref_atd(1:3);
rSun_aux_i = ref_atd(10:12);
rNadir_i = -eX_orb(1:3)'/Mag(eX_orb(1:3)');
v_i = eX_orb(4:6)'/Mag(eX_orb(4:6)');

p_fss = scdata.p_fss;
p_uhf = scdata.p_uhf;
p_aux = scdata.p_aux;
p_drag = scdata.p_drag;
p_uhf_aux = scdata.p_uhf_aux;

q = eX_atd(4:7);

if flag_AC == 1
    
    if flag_fss == 1
        
        [q_err,angle_err] = Cal_AtdErr_2(rSun_i,rSun_aux_i,p_fss,p_aux,q);
        w_err = eX_atd(1:3)';
        
    else
        
        [q_err,angle_err] = Cal_AtdErr_1(rSun_i,p_fss,q);
        w_err = eX_atd(1:3)';
        
    end
    
    if q_err(1) < 0
        q_err = -q_err;
    end
    
    [ T ] = Ctrl_Atd_PD( q_err,w_err,T_max,ctrlgain_atd_PD(:,1) );
    %         angle_err
    
elseif flag_AC == 2
    
    %         [q_err,angle_err] = Cal_AtdErr_1(rNadir_i,p_uhf,q);
    [q_err,angle_err] = Cal_AtdErr_2(rNadir_i,v_i,p_uhf,p_uhf_aux,q);
    w_err = eX_atd(1:3)';
    
    if q_err(1) < 0
        q_err = -q_err;
    end
    
    [ T ] = Ctrl_Atd_PD( q_err,w_err,T_max,ctrlgain_atd_PD(:,2) );
    %         angle_err
    
elseif flag_AC == 3
    
    %         if flag_fss == 1
    %
    %             v_aux_i = cross(rSun_i,v_i);
    %             p_drag_aux = cross(p_fss,p_drag);
    %             [q_err,angle_err] = Cal_AtdErr_2(rSun_i,v_aux_i,p_fss,p_drag_aux,q);
    %             w_err = eX_atd(1:3)';
    %
    %         else
    %
    %             [q_err,angle_err] = Cal_AtdErr_1(v_i,p_drag,q);
    %             w_err = eX_atd(1:3)' - 0.5*pi/180*p_drag;
    %
    %         end
    %
    %         if q_err(1) < 0
    %             q_err = -q_err;
    %         end
    %
    %         [ T ] = Ctrl_Atd_PD( q_err,w_err,T_max,ctrlgain_atd_PD );
    % %         angle_err
    
    [q_err,angle_err] = Cal_AtdErr_2(rNadir_i,v_i,-p_uhf,p_uhf_aux,q);
    w_err = eX_atd(1:3)';
    
    if q_err(1) < 0
        q_err = -q_err;
    end
    
    [ T ] = Ctrl_Atd_PD( q_err,w_err,T_max,ctrlgain_atd_PD(:,2) );
    
elseif flag_AC == 4
    q_d = quat_targ;
    q_err = quatnormalize(quatmultiply([q_d(1,1) -q_d(1,2:4)],q));
    %         q_err = quatnormalize(quatmultiply(q,[q_d(1,1) -q_d(1,2:4)]));
    %         q_err = quatnormalize(quatmultiply(q_d,q_inv(q)'));
    w_err = eX_atd(1:3)';
    
    if q_err(1) < 0
        q_err = -q_err;
    end
    
    [ T ] = Ctrl_Atd_PD( q_err,w_err,T_max,ctrlgain_atd_PD(:,submode+2) );
    %         [ T ] = Ctrl_Atd_PD( w_err,q_inv(q_err)',T_max,ctrlgain_atd_PD(:,submode+2) );
    
else
    
    T = zeros(3,1);
    angle_err = 0;
    q_err = zeros(1,4);
    
end
end