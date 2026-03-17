function dot = HW3_stateeq(t, state, I)
omx = state(1);
omy = state(2);
omz = state(3);
psi = state(4);
theta = state(5);
phi = state(6);

dot(1, 1) = (I(2) - I(3))/I(1)*omy*omz; 
dot(2, 1) = (I(3) - I(1))/I(2)*omz*omx;
dot(3, 1) = (I(1) - I(2))/I(3)*omx*omy;
dot(4, 1) = (omy*sin(phi) + omz*cos(phi))*sec(theta);
dot(5, 1) = omy*cos(phi) - omz*sin(phi);
dot(6, 1) = omx + omy*sin(phi)*tan(theta) + omz*cos(phi)*tan(theta);