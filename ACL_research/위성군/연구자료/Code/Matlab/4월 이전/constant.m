clear all
clc

%% 기본상수
r = 6378.14               % km
mu = 398600.44            % km^3/s^2
j_2 = 0.00108263
al = 570
a = al + r                 % km
i = deg2rad(102.9)         % rad
e = 0
p = a*(1-e^2)
w_e = 0.0000729            % rad/s 
n = sqrt(mu/a^3)           % rad/s
G = 6.674 * 10^(-11)
M = 5.972*10^24

M0_dot = -3*n*r^2*j_2*sqrt(1-e^2)*(3*sin(i)^2-2)/(4*p^2);       % rad/s
M_dot = M0_dot + n
w_dot = (3/2)*(j_2)*(r/p)^2*n*(2-(5/2)*(sin(i)^2));             % rad/s
Omega_dot = -(3/2)*(j_2)*(r/p)^2*n*(cos(i)) 

TS = 2*pi/(w_dot + M0_dot);       % s 
TG = 2*pi/(w_e - w_dot);          % s 

T = 2*pi/n         % s



ae0 = 86400/T;  % 섭동 고려 못할 때 (8058.98)

ae1 = (w_dot + M_dot) / (w_e - Omega_dot);  %tau = TG/TS %(NP/ND) 섭동 고려 (8056.11)


%d_Omega = 2*pi / (365.25*86400) 

%% 지구 상수

% w_dot = (3/2)*(j_2)*(r/p)^2*n*(2-(5/2)*(sin(i)^2))    % rad/s
% Omega_dot = -(3/2)*(j_2)*(r/p)^2*n*(cos(i))    % rad/s
% M_dot = n*(1-(3/2)*(j_2)*(r/p)^2)*sqrt(1-e^2)*((3/2)*(sin(i))^(2)-1)    % rad/s
% 
% TS = 2*pi/(w_dot + M_dot)    % s 
% TG = 2*pi/(w_e - Omega_dot)    % s 
% 
% 
% tau = TG/TS %(NP/ND)



% omega_k = (n_k*(2*pi)*N_D) / L + omega_0


% %% SSO 궤도 정보(궤도 경사각)
% cosi = -0.2233
% sini = 0.9748
% inc_deg = acos(cosi) * 180 / pi

