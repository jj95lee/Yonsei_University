function K_dot = K_eq(t, K)
K_dot(1) = 2*(K(2)^3 - 2*K(2) - 1);
K_dot(2) = 2*K(2)*K(3) - K(1) + K(2) - 2*K(3);
K_dot(3) = 2*K(3)^2 - 2*K(2) + 2*K(3) - 1;