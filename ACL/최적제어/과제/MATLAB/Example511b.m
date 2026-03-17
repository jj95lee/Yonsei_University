% Initial/final times and boundary conditions
t0 = 0; tf = 2;
x0 = [0; 0];


% A matrix
A = [0 1 0 0; 0 -1 0 -1; 0 0 0 0; 0 0 -1 1];
expA_f = expm(A * tf);
Lambda0 = [expA_f(1, 3)+5*expA_f(2, 3) expA_f(1, 4)+5*expA_f(2, 4); ...
           5*expA_f(3, 3)-expA_f(4, 3) 5*expA_f(3, 4)-expA_f(4, 4)]^-1 ...
           * [15; 0];


% t, x, Lambda, u matrices
t = [t0:(tf - t0)/100:tf]';
x1 = zeros(size(t));
x2 = zeros(size(t));
Lambda1 = zeros(size(t));
Lambda2 = zeros(size(t));
u = zeros(size(t));

for i = 1:length(t)
    xL_t = expm(A * t(i)) * [x0; Lambda0];
    x1(i) = xL_t(1);
    x2(i) = xL_t(2);
    Lambda1(i) = xL_t(3);
    Lambda2(i) = xL_t(4);
end
u = -Lambda2; 

plot(t, x1, t, x2, t, Lambda1, t, Lambda2, t, u)
title('Example 5.1-1c')
xlabel('Time'); ylabel('State, costate and control')
legend({'x1', 'x2', '\Lambda1', '\Lambda2', 'u'});
grid on