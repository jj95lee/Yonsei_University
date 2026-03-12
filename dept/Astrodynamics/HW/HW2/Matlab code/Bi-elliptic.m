% In this example we compare be-elliptic and Hohmann transfer.
% The total speed change that required for spacecraft transfer from
% geocentric circular orbit with radius Ri to a higher altitude Rf
clc;
clear all;
Rf  =  42159;   % [km] Final circular orbit
Ri  =  6569;     % [km] Initial circular orbit
Rb  =  54200;   % [km] Apogee of the transfer ellipse
mu  =  398600;   % [km^3/s^2] Earth’s gravitational parameter

% For initial circular orbit
Vc = (mu/Ri)^0.5;
a   = Rf/Ri;
b   = Rb/Ri;

dV_BE = Vc*((2*(a+b)/(a*b))^0.5 - (1+1/a^0.5) - ((2/(b +b^2))^0.5*(1-b)));
% Semimajor axes of the first transfer ellipse
a1 = (Ri + Rb)/2;
% Semimajor axes of the second transfer ellipse
a2 = (Rf + Rb)/2;
t_BE = pi/(mu)^0.5*(a1^1.5+a2^1.5);     %[s]
fprintf('Total speed change = %4.4f [km/s]\n',dV_BE);
fprintf('Time required for transfer = %4.2f [hours]\n',t_BE/3600);