r_int=[5328.7862 4436.1273 101.4720];
v_int=[-4.864779 5.816486 0.240163];

r_tgt=[6697.4756 1794.5831 0];
v_tgt=[-1.962372 7.323674 0];
t=6000;
mu=398600;
[p,a,e,i,w,argp,nu,m,arglat,truelon,longper] = rv2coe (r_tgt,v_tgt,mu);