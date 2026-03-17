% Initial/final times and boundary conditions
t0 = 0; tf = 1;
x0 = 1; xf = 0;


% A matrix
A = [0 -1; -1 0];
expA_f = expm(A * tf);
Lambda0 = expA_f(1, 2)^-1 * (xf - expA_f(1, 1) * x0);


% t, x, Lambda, u matrices
t = [t0:(tf - t0)/100:tf]';
x = zeros(size(t));
Lambda = zeros(size(t));
u = zeros(size(t));

for i = 1:length(t)
    xL_t = expm(A * t(i)) * [x0; Lambda0];
    x(i) = xL_t(1);
    Lambda(i) = xL_t(2);
end
u = -Lambda; 

plot(t, x, t, u)
title('Example on p. 58')
xlabel('Time'); ylabel('State and control')
legend({'x', 'u'});
grid on