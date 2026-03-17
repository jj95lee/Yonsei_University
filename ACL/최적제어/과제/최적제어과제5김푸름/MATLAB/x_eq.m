function x_dot = x_eq(t, x)

global a K_history

x_dot = (a - 2 * K_history(round(t * 100 + 1))) * x;
