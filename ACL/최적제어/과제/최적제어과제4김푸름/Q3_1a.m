A = [0 1 0 0; -1 -1 0 -1; -1 0 0 1; 0 -1 -1 1];
eA1 = expm(A * 1);
L0 = eA1(1:2, 3:4)^-1 * ([0; 0] - eA1(1:2, 1:2) * [1; 0]);

t = 0:0.01:1;
t = t';
x1 = zeros(size(t));
x2 = zeros(size(t));
L1 = zeros(size(t));
L2 = zeros(size(t));
u = zeros(size(t));

for i = 1:length(t)
    Xt =  expm(A * t(i)) * [1; 0; L0];
    x1(i) = Xt(1);
    x2(i) = Xt(2);
    L1(i) = Xt(3);
    L2(i) = Xt(4);
end
u = -L2; 

plot(t, x1, t, x2, t, u)
title('Problem 3, Case 1a)')
xlabel('Time'); ylabel('State and control')
legend({'x_{1}', 'x_{2}', 'u'});
grid on