% function plot_FSS(option)
    fov=120;
%     switch option
%         case 1 %alpha data
            for i=1:360
                alpha=i/180*pi;
                e=0/180*pi;
                [alpha_cal(i), e_cal_dummy]=FSS(alpha,e,fov);
                alpha_cal_error(i)=alpha_cal(i)-i;
            end
%         case 2
            for i=fov/2:90
                alpha=10/180*pi;
                e=i/180*pi;
                [alpha_cal_dummy, e_cal(i)]=FSS(alpha,e,fov);
                e_cal_error(i)=e(i)-i;
            end
%     end
    
    figure()
    plot(1:360,alpha_cal,'linewidth',1.6);
    xlabel('true (deg)');
    ylabel('modeled alpha (deg)')
    
    figure()
    plot(fov/2:90,e,'linewidth',1.6);
    xlabel('true (deg)');
    ylabel('modeled epsilon (deg)')
    
    figure()
    plot(1:360,alpha_cal_error,'linewidth',1.6);
    xlabel('true (deg)');
    ylabel('modeled alpha error (deg)')
    
    figure()
    plot(fov/2:90,e_cal_error,'linewidth',1.6);
    xlabel('true (deg)');
    ylabel('modeled epsilon error (deg)')
    

    
% end
    