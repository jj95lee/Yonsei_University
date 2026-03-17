x(1) = 1;
x(2) = x(1) - f1(x(1))/f1dot(x(1));
idx = 2;
while abs((x(idx) - x(idx - 1)) / x(idx)) > 0.0001
    x(idx + 1) = x(idx) - f1(x(idx))/f1dot(x(idx));
    idx = idx + 1;
end