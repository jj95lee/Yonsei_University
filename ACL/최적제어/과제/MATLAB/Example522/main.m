global A B Q R
A = [0 1; 2 -1];
B = [0; 1];
Q = [2 0; 0 1];
R = 0.5;
x0 = [-4; 4];
t = [0:0.1:15]';

syms k11 k12 k22
K = [k11 k12; k12 k22];
K_dot = -Q - A' * K - K * A + K * B * (1/R) * B' * K;

K_history = ode45('K_eq', t, [0 0 0]);
x_history = ode45('x_eq', t, x0);
u_history = -(1/R) * B' * (K_history(:, 2) .* x_history(:, 1) + K_history(:, 3) .* x_history(:, 2));