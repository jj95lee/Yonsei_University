r_int=[5328.7862 4436.1273 101.4720];
v_int=[-4.864779 5.816486 0.240163];

r_tgt=[6697.4756 1794.5831 0];
v_tgt=[-1.962372 7.323674 0];
t=6000;
mu=398600;
[p,a,e,i,w,argp,nu,m,arglat,truelon,longper ] = rv2coe (r_tgt,v_tgt, mu);


[r_TGT,v_TGT,errork] = kepler(r_tgt,v_tgt,t);

[p,a,e,i,w,argp,nu1,m,arglat,truelon,lonper ] = rv2coe (r_TGT,v_TGT, mu);


m=mu*t^2/((2*sqrt(norm(r_int)*norm(r_TGT))*cos((nu1-nu)/2))^3);
l=(norm(r_int)+norm(r_TGT))/(4*sqrt(norm(r_int)*norm(r_TGT))...
*cos((nu1-nu)/2))-1/2;

q0=2*atan(sqrt((1-norm(e))/1+norm(e))*tan(nu/2));
q=2*atan(sqrt((1-norm(e))/1+norm(e))*tan(nu1/2));

X=(q-q0-sin(q-q0))/sin((q-q0)/2)^3;
Y=[1, -1, 0, -m*X];
Y=roots(Y);
y=Y(3);

p=y^2*norm(r_int)^2*norm(r_TGT)^2*sin(nu1-nu)^2/(mu*t^2);

f=1-r_TGT/p*(1-cos(nu1-nu));
g=norm(r_int)*norm(r_TGT)*sin(nu1-nu)/(sqrt(mu*p));
ff=sqrt(1/p)*tan((nu1-nu)/2)*((1-cos(nu1-nu))/p-1/norm(r_int)-1/norm(r_TGT));
gg=1-r_int/p*(1-cos(nu1-nu));

v_0=(r_TGT-f.*r_int)./g;
v=(gg.*r_TGT-r_int)./g;

delv_a=v_0-v_int;
delv_b=v_tgt-v;


if dot(r_int,v_0)<0
    if dot(r_tgt,v)>0
    mu=398600;
    a=-mu/2*(norm(v_0)^2/2-mu/norm(r_int));
    h=cross(r_int,v_0);
    p=norm(h)^2/mu;
    e=sqrt((a-p)/a);
    r_p=a*(1-e);
    if r_p<=6378.136
        fprintf('Escape Earth\n')
    else
        fprintf('We are safe\n')
    end
    else
        fprintf('We are safe2\n')
    end
else
fprintf('We are safe1\n')
end



