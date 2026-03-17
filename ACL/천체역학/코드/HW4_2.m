om0 = (2*pi)/86400;
Ix = 1200
Iy = 1000;
Iz = 800;

A = [0 1; 3*om0^2*(Iz - Ix)/Iy 0]

eig(A)

%%
B = [0              om0                 0               1;
     -om0           0                   1               0;
     0              3*om0^2*(Iz-Iy)/Ix  0               om0*(Iy-Iz)/Ix;
     0              0                   om0*(Ix-Iy)/Iz  0];

eig(B)

%%
omx = om0 * cos(0.1) * sin(0)
omy = om0 * (cos(0.1) * cos(0) + sin(0.1) * sin(0.1) * sin(0))
omz = om0 * (-sin(0.1) * cos(0) + cos(0.1) * sin(0.1) * sin(0))


%%
func = @(t, x) ode_att(t, x, [Ix Iy Iz], om0);
[list_t, list_x] = ode45(func, 0:86400, [omx omy omz 0 0.1 0.1]);

f1 = figure(1);
set(f1, 'Position', [100 100 900 900])

subplot(3, 2, 1)
plot(0:24/86400:24, list_x(:, 1))
grid on

xlabel('Time (h)')
xlim([0 24])
set(gca, 'xtick', [0:24:240])
ylabel('\omega_x (rad/sec)')


subplot(3, 2, 3)
plot(0:24/86400:24, list_x(:, 2))
grid on

xlabel('Time (h)')
xlim([0 24])
set(gca, 'xtick', [0:24:240])
ylabel('\omega_y (rad/sec)')


subplot(3, 2, 5)
plot(0:24/86400:24, list_x(:, 3))
grid on

xlabel('Time (h)')
xlim([0 24])
set(gca, 'xtick', [0:24:240])
ylabel('\omega_z (rad/sec)')


subplot(3, 2, 2)
plot(0:24/86400:24, list_x(:, 4))
grid on

xlabel('Time (h)')
xlim([0 24])
set(gca, 'xtick', [0:24:240])
ylabel('\psi (rad)')


subplot(3, 2, 4)
plot(0:24/86400:24, list_x(:, 5))
grid on

xlabel('Time (h)')
xlim([0 24])
set(gca, 'xtick', [0:24:240])
ylabel('\theta (rad)')


subplot(3, 2, 6)
plot(0:24/86400:24, list_x(:, 6))
grid on

xlabel('Time (h)')
xlim([0 24])
set(gca, 'xtick', [0:24:240])
ylabel('\phi (rad)')

%%
func = @(t, x) ode_att(t, x, [Ix Iy Iz], om0);
[list_t, list_x] = ode45(func, 0:864000, [omx omy omz 0 0.1 0.1]);

f2 = figure(2);
set(f2, 'Position', [100 100 900 900])

subplot(3, 2, 1)
plot(0:24/86400:240, list_x(:, 1))
grid on

xlabel('Time (h)')
xlim([0 240])
set(gca, 'xtick', [0:24:240])
ylabel('\omega_x (rad/sec)')


subplot(3, 2, 3)
plot(0:24/86400:240, list_x(:, 2))
grid on

xlabel('Time (h)')
xlim([0 240])
set(gca, 'xtick', [0:24:240])
ylabel('\omega_y (rad/sec)')


subplot(3, 2, 5)
plot(0:24/86400:240, list_x(:, 3))
grid on

xlabel('Time (h)')
xlim([0 240])
set(gca, 'xtick', [0:24:240])
ylabel('\omega_z (rad/sec)')


subplot(3, 2, 2)
plot(0:24/86400:240, list_x(:, 4))
grid on

xlabel('Time (h)')
xlim([0 240])
set(gca, 'xtick', [0:24:240])
ylabel('\psi (rad)')


subplot(3, 2, 4)
plot(0:24/86400:240, list_x(:, 5))
grid on

xlabel('Time (h)')
xlim([0 240])
set(gca, 'xtick', [0:24:240])
ylabel('\theta (rad)')


subplot(3, 2, 6)
plot(0:24/86400:240, list_x(:, 6))
grid on

xlabel('Time (h)')
xlim([0 240])
set(gca, 'xtick', [0:24:240])
ylabel('\phi (rad)')