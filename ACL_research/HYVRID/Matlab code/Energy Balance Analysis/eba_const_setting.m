%% Gravitation
mu = 398600; % Earth Gravitational Standard Constant, km3/s2
J2 = 0.00108262; % J2 Perturbation Coefficient, -
Re = 6378.1368; % Radius of Earth, km

%% Orbit and Attitude
H = 500; % Altitude, km
Tepoch = [2025, 1, 1]; % Launch Date, year/month/day
MLTAN = [27, 30, 0]; % Local Time on Ascending Node, hour/minute/sec
Flag_att0 = 1; % Nominal Attitude Orientation, 1 = Nadir-pointing, 2 = Sun-pointing

%% PV Cell
Npvs = 20;
Npvp = 1;
Vmp_eol = 2.32;
Imp_eol = 0.490; 

%% Battery
SOC0 = 70; % Initial Battery SOC, %
Bat_cap = 1.241; % Battery Capcity per Cell, Ah
% OCV: Open Circuit Voltage, V
OCV = [3
    3.61
    3.67
    3.7
    3.72
    3.74
    3.76
    3.77
    3.78
    3.79
    3.8
    3.82
    3.85
    3.88
    3.92
    3.95
    3.98
    4.03
    4.07
    4.14
    4.2
    ];
Ri_cell = (140+190)/2*1e-3; % Internal Impedence of Battery, mOhm
Ns = 2; % Number of Battery Cells in Series Connection, -
Np = 4; % Number of Battery Cells in Parallel Connection, -
V_ulim = 8.26; % Upper Limit to Prevent Battery Over Charging, V

%% Etc.
eff_tmp = 90; % Efficiency for Charging due to Temperature of PV Cell, % 
eff_id = 90; % Efficiency for Charging due to Inherent Degradation, %
eff_bcr = 80; % Conversion Efficiency for Charging Battery, %
eff_pcm = 90; % Conversion Efficiency for Regulation, %