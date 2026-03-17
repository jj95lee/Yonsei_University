function Phi = STM_2body_inplane(omega_0, t)

o = omega_0;

Phi = [4 - 3*cos(o*t)       sin(o*t)/o              0       2*(1-cos(o*t))/o;
       3*o*sin(o*t)         cos(o*t)                0       2*sin(o*t);
       6*sin(o*t)-6*(o*t)   2*(-1 + cos(o*t))/o     1       4*sin(o*t)/o - 3*t;
       6*o*(-1 + cos(o*t))  -2*sin(o*t)             0       -3 + 4*cos(o*t)];