function x_dot = x_eq(t, x)

global a T L

x_dot = (a - 2 * K(a, T, L, t)) * x;
