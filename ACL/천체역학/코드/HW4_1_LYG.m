om0 = 0.0011;
Ix = 1400;
Iy = 1700;
Iz = 1000;

A = [0 1; 3*om0^2*(Iz - Ix)/Iy 0]

eig(A)

%%
B = [0              om0                 0               1;
     -om0           0                   1               0;
     0              3*om0^2*(Iz-Iy)/Ix  0               om0*(Iy-Iz)/Ix;
     0              0                   om0*(Ix-Iy)/Iz  0];

eig(B)

%%
omx = om0 * cos(0.1) * sin(0.1)
omy = om0 * (cos(0.1) * cos(0.1) + sin(0.1) * sin(0.1) * sin(0.1))
omz = om0 * (-sin(0.1) * cos(0.1) + cos(0.1) * sin(0.1) * sin(0.1))


%%
func = @(t, x) ode_att(t, x, [Ix Iy Iz], om0);
[list_t, list_x] = ode45(func, 0:86400, [omx omy omz 0.1 0.1 0.1]);

f1 = figure(1);
set(f1, 'Position', [100 100 900 900])

subplot(3, 2, 1)
plot(0:24/86400:24, list_x(:, 1))
grid on

xlabel('t(h)')
xlim([0 12])
ylim([-0.0005 0.0005])
set(gca, 'xtick', [0:6:24])
ylabel('Wx(rad/s)')


subplot(3, 2, 3)
plot(0:24/86400:24, list_x(:, 2))
grid on

xlabel('t(h)')
xlim([0 12])
ylim([0.0009 0.0013])
set(gca, 'xtick', [0:6:24])
ylabel('Wy(rad/s)')


subplot(3, 2, 5)
plot(0:24/86400:24, list_x(:, 3))
grid on

xlabel('t(h)')
xlim([0 12])
ylim([-0.0003 0.0003])
set(gca, 'xtick', [0:6:24])
ylabel('Wz(rad/s)')


subplot(3, 2, 2)
plot(0:24/86400:24, list_x(:, 4))
grid on

xlabel('t(h)')
xlim([0 12])
ylim([-0.4 0.4])
set(gca, 'xtick', [0:6:24])
ylabel('\psi(rad)')


subplot(3, 2, 4)
plot(0:24/86400:24, list_x(:, 5))
grid on

xlabel('t(h)')
xlim([0 12])
ylim([-0.2 0.2])
set(gca, 'xtick', [0:6:24])
ylabel('\theta(rad)')


subplot(3, 2, 6)
plot(0:24/86400:24, list_x(:, 6))
grid on

xlabel('t(h)')
xlim([0 12])
ylim([-0.2 0.2])
set(gca, 'xtick', [0:6:24])
ylabel('\phi(rad)')