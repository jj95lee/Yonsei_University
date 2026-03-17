function Pmp = cal_TJSC3G30_mp(G, Ns, TC)
% Paramter
Gref = 1367; % Solar Radiation, W/m2
q = 1.60217646e-19; % Electron Charge, C = A.s
k = 1.3806503e-23; % Boltzmann Constant, J/K = kg.m2/s2/K

A = 30.18*1e-4; % Cell Area, m2

% Temperature Gradient Coefficient
Kv = -6.2/2700/1000*Ns; % V/deg.C
Ki = 0.36/520.2/1000; % A/deg.C

Voc = 2700/1000*Ns; % V
Isc = 520.2/1000; % A

% Pindado, 2016
Rs = 7.95*1e-2; % Equivalent Series Resistance, Ohm
Rp = 2.62*1e3; % Shunt Resistance, Oh
alpha = 3.72; % Diode Ideality Constant

% PV Cell Current Calculation
TK = TC + 273.15; % Kelvin

Vt = Ns*k*TK/q; % Thermal Voltage

VocT = Voc + Voc*Kv*(TC - 28);
IscT = Isc + Isc*Ki*(TC - 28);

I0 = IscT/(exp(VocT/alpha/Vt) - 1);
Ipv = IscT*(G/Gref);

Iini = Isc;

V = 0:Voc/100:Voc;
for k = 1:length(V)
    Ir(k) = Ipv - I0*(exp((V(k) + Rs*Iini)/alpha/Vt) - 1) - (V(k) + Rs*Iini)/Rp;
    P(k) = Ir(k)*V(k);
end

Pmp = max(P);