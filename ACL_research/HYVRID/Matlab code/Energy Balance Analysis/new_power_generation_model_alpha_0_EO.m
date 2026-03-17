close all
clear all
clc

%% Call EBA Constant Setting
eba_const_setting

%% Orbit Dynamics
sma = Re + H; % Semi-major Axis, km
ecc = 0; % Eccentricity
anr = 360/365.25*pi/180/86400; % Ascending Node Angular Rate, deg/sec

Tvernal = [2025, 3, 20]; % Date of Vernal Equinox (* Set Same Year with Epoch Time)
if Tepoch(1) > Tvernal(1)
    Tepoch(1) = Tepoch(1) - floor(Tepoch(1) - Tvernal(1));
    DayDiff = char(between(datetime(Tepoch), datetime(Tvernal), 'Days')); 
else
    DayDiff = char(between(datetime(Tepoch), datetime(Tvernal), 'Days'));
end
SunLong = str2num(DayDiff(2:(strlength(DayDiff)-1)))*360/365.25; % Ecliptic True Solar Longitude, deg
eps = 23.45; % Obliquity of the Ecliptic, deg 
UTC = [MLTAN(1)-9, 30, 0]; % Coordinated Universal Time of Ascending Node, hour/minute/sec 
mltan; % Run Convertor, MLTAN to RAAN
RAAN = rtd*raan; % *Calculate RAAN using mltan.m, deg
incl = acosd(-anr/(3/2*J2*(Re/sma/(1 - ecc^2))^2*(mu/sma^3)^0.5)); % Inclination, deg
beta = asind(cosd(SunLong)*sind(RAAN)*sind(incl) - sind(SunLong)*cosd(eps)*cosd(RAAN)*sind(incl) + sind(SunLong)*sind(eps)*cosd(incl)); % Beta Angle, deg

%% Initial Value
ang_in(1) = 270 - acosd(Re/(Re+H));
ang(1) = 0;
t(1) = 0;
j = 1;
alpha = acosd(Re/(Re+H))*0;

%% Position in Orbit as Angle

while sqrt(mu/sma)/sma*180/pi*t(j) <= 360
    ang_in(j, 1) = sqrt(mu/sma)/sma*180/pi*t(j) - alpha + 270;
    ang(j,1) = sqrt(mu/sma)/sma*180/pi*t(j);
    if   ang_in(j, 1) > 479.4
       ang_in(j, 1) = 0;
    elseif ang_in(j, 1) < 292
        ang_in(j, 1) = 0;
    elseif ang_in(j, 1) > 360
        ang_in(j, 1) = (sqrt(mu/sma)/sma*180/pi*t(j) - alpha) - 90;
    end

    t(j + 1) = t(j) + 1;
    j = j + 1;
end

[r, c] = size(ang_in(:, 1));

%% Calculate Maximum Power Generation in Each Cell (or Panel)
for i = 1:r
   
        if ang_in(i) > 0
            Pmp(i, 1) = 1.135*(6*3+6*3+4*3+4*3)/(3+3+3+3+1+1);           
        else
            Pmp(i, 1) = 0;
        end
    
    Pgen(i, :) = Pmp(i, 1);
    
end

%% Figure
figure()
set(gcf,'Color','w')
plot(Pgen(:, 1)); hold on
% plot(Pmp(:, 2)); hold on
% plot(Pmp(:, 3)); hold on
% plot(Pmp(:, 4)); hold on
% plot(Pmp(:, 5)); hold on
% plot(Pmp(:, 6))
grid on
% ylim([0 100])
xlabel('Time (sec)')
ylabel('Maximum Power Generation (W)')
% legend('P_{Total}', 'P_{+X}', 'P_{-X}', 'P_{+Y}', 'P_{-Y}', 'P_{+Z}', 'P_{-Z}')

%% Set Paramters to Estimate SOC
t = [0:1:size(Pgen)-1]';
t_end = max(t);
fc = 0.01;
SOC_init = SOC0/100;
SOC = [0:5:100];
