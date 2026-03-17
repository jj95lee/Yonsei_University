mu_E = 3.986004418E14;
R_580 = (6.3781 + 0.5800) * 10^6;
R_600 = (6.3781 + 0.6000) * 10^6;
t_end = 1800;



state_original_0 = [R_580, 0, 0, 0, R_580 * sqrt(mu_E / R_580^3), 0];
state_deputy_0 = [R_600, 0, 0, 0, R_600 * sqrt(mu_E / R_600^3), 0];
state_original_0_pert = [R_580, 0, 0, 0 + 23.054, R_580 * sqrt(mu_E / R_580^3) + 0.088, 0];


opt_ode45 = odeset('RelTol', 1e-12);
tB = @(t, X)twoBody(t, X);


[~, state_original] = ode45(tB, [0 t_end], state_original_0, opt_ode45);
[~, state_deputy] = ode45(tB, [0 t_end], state_deputy_0, opt_ode45);
[~, state_original_pert] = ode45(tB, [0 t_end], state_original_0_pert, opt_ode45);

hold on
plot(state_deputy(:, 1), state_deputy(:, 2));
plot(state_original(:, 1), state_original(:, 2));
plot(state_original_pert(:, 1), state_original_pert(:, 2));

function Xdot = twoBody(t, X)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% t: Time since integration (s)
% X: J2000 heliocentric ecliptic state (X, Y, Z, VX, VY, VZ) (m, m, m, m s^-1, m s^-1, m s^-1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mu_E = 3.986004418E14; % Sun standard gravitational parameter

r = sqrt(X(1)^2 + X(2)^2 + X(3)^2);
Xdot = [X(4); X(5); X(6);
        (-mu_E / r^3) * X(1);
        (-mu_E / r^3) * X(2);
        (-mu_E / r^3) * X(3)];
end