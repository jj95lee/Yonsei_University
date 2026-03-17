function xdot = ode_att(t, x, I, om0)

xdot(1, 1) = (I(2) - I(3)) * x(2)*x(3) - 3*om0^2*(I(2) - I(3))*sin(x(6))*cos(x(6))*cos(x(5))^2;
xdot(1, 1) = xdot(1) / I(1);

xdot(2, 1) = (I(3) - I(1)) * x(3)*x(1) + 3*om0^2*(I(3) - I(1))*cos(x(6))*sin(x(5))*cos(x(5));
xdot(2, 1) = xdot(2) / I(2);

xdot(3, 1) = (I(1) - I(2)) * x(1)*x(2) + 3*om0^2*(I(1) - I(2))*sin(x(6))*sin(x(5))*cos(x(5));
xdot(3, 1) = xdot(3) / I(3);

wxr = x(1) - om0*cos(x(5))*sin(x(4));
wyr = x(2) - om0*(cos(x(6)*cos(x(4)) + sin(x(6))*sin(x(5))*sin(x(4))));
wzr = x(3) - om0*(-sin(x(6))*cos(x(4)) + cos(x(6))*sin(x(5))*sin(x(4)));

xdot(4, 1) = (wyr*sin(x(6)) + wzr*cos(x(6)))*sec(x(5));
xdot(5, 1) = wyr*cos(x(6)) - wzr*sin(x(6));
xdot(6, 1) = wxr + wyr*sin(x(6))*tan(x(5)) + wzr*cos(x(6))*tan(x(5));

