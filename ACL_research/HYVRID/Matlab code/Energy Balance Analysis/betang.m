function [bet] = betang(beta,theta, ANG, camorant,flag_Pin )
    rot1 = [     1 ,            0,               0    ;           
                 0 ,         cosd(ANG) ,    sind(ANG) ;
                 0 ,        -sind(ANG) ,    cosd(ANG)       ] ;
    rot2 = [ cosd(beta) ,         0 ,       -sind(beta) ;
                 0  ,            1 ,             0 ;
             sind(beta) ,         0 ,        cosd(beta)       ] ;
    rot3 = [ cosd(ANG) ,sind(ANG) ,        0 ;
            -sind(ANG) ,cosd(ANG) ,        0 ;
                 0 ,            0,               1          ] ;
    rot4 = [ cosd(90) ,sind(90) ,        0 ;
            -sind(90) ,cosd(90) ,        0 ;
                 0 ,            0,               1          ] ;
    sun   = [ cosd(beta) ,       0,            sind(beta)   ] ;
    x(:,1)  = [-cosd(theta); 0; sind(theta)]; %DSA1 3
    x(:,2)  = [ cosd(theta); 0;-sind(theta)]; %DSA2 3
    x(:,3)  = [ cosd(theta); 0; sind(theta)]; %DSA3 3
    x(:,4)  = [-cosd(theta); 0;-sind(theta)]; %DSA4 3
    x(:,5)  = [  0;  1 ; 0] ; %BSAZM 1
    x(:,6) = [  0;  0 ;-1] ; %BSAYM 2
    x(:,7) = [  0;  0 ; 1] ; %BSAYP 2
    x(:,8) = [  -1; 0 ; 0] ; %BSAXP 2
    x(:,9) = [   1; 0 ; 0] ; %BSAXM 2
    if flag_Pin == 1
        if camorant == 1
            for i = 1:9
                x(:,i) = rot3 * x(:,i) ;
                c(i) = dot(x(:,i),sun);
                bet(i,1) = acosd(c(i));
            end
        else
            for i = 1:9
                x(:,i) = rot3 *rot4 * x(:,i) ;
                c(i) = dot(x(:,i),sun);
                bet(i,1) = acosd(c(i));
            end
        end
    else  
        for i = 1:9
            x(:,i) = rot2 * x(:,i) ;
            c(i) = dot(x(:,i),sun);
            bet(i,1) = acosd(c(i)); 
        end
end