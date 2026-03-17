global a T L x0 K_history
for L = [5, 0.5, 0.05]
    a = -0.2;
    T = 15;
    x0 = 5;
    t = [0:0.01:T]';

    c_1 = 2 * T + log(L / (L - a)) / a;
    K_history = a * exp(a * c_1)./(exp(a * c_1) - exp(2 * a * t));

    [tt, x_history] = ode45('x_eq', t, x0);
    u_history = - 2 * K_history .* x_history;

    figure
    plot(t, K_history, t, u_history, t, x_history)
    title(['L = H = ' num2str(L)]);
    xlabel('Time (t)');
    legend('K(t)', 'u(t)', 'x(t)');
end