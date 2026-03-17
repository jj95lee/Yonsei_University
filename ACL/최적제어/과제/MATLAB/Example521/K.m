function Kresult = K(a, T, L, t)

temp1 = 1./( exp(-a * (T - t)) - (L/a)*(exp(-a * (T - t)) - exp(a * (T - t))) );
Kresult = temp1 .* ( L * exp(a*(T-t)) );