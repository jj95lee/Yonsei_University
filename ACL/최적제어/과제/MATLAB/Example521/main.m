global a T L x0
a = -0.2;
T = 15;
L = 5;
x0 = 5;
t = [0:0.01:T]';

K_history = K(a, T, L, t);

[tt, x_history] = ode45('x_eq', t, x0);
u_history = - 2 * K_history .* x_history;

plot(t, K_history, t, u_history, t, x_history)
xlabel('Time (t)');
legend('K(t)', 'u(t)', 'x(t)');