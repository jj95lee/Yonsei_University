clc
clear

global mu req hn1 hn2 hn3 


req = 6378;    % Radius of earth
r1 = 6670;         %Initial radii (km)
r2 = 42164;        %Final radii (km)
a1 = 8893;

% compute "normalized" radii

hn1 = sqrt(2.0 * r2 / (r2 + r1));
hn2 = sqrt(r1 / r2);
hn3 = sqrt(2.0 * r1 / (r2 + r1));

% compute "local circular velocity" of initial and final orbits (km/sec)
mu = 398600.436233;
v1 = sqrt((2mu / r1)-(mu / a1);
v2 = sqrt(mu / r2);

% compute transfer orbit semimajor axis (kilometers)
smat = 0.5 * (r1 + r2);

% compute transfer orbit eccentricity (non-dimensional)
ecct = (max(r1, r2) - min(r1, r2)) / (r1 + r2);

% compute transfer orbit perigee and apogee radii and velocities
rp = smat * (1.0 - ecct);
ra = smat * (1.0 + ecct);
vt1 = sqrt(2.0 * mu * ra / (rp * (rp + ra)));
vt2 = sqrt(2.0 * mu * rp / (ra * (rp + ra)));

% compute transfer orbit period (seconds)

taut = 2.0 * pi * sqrt(smat^3 / mu);
tof = 0.5 * taut;

dv1 = vt1 - v1;
dv2 = v2 - vt2;


% print results
clc; home;

fprintf('\nHohmann Orbit Transfer Analysis');
fprintf('\n-------------------------------\n\n');

fprintf('initial orbit radius              %10.4f kilometers \n\n', r1);

fprintf('initial orbit velocity            %10.4f meters/second \n\n\n', 1000.0 * v1);

fprintf('final orbit radius                %10.4f kilometers \n\n', r2);

fprintf('final orbit velocity              %10.4f meters/second \n', 1000.0 * v2);

fprintf('first delta-v                     %10.4f meters/second \n\n', 1000.0 * dv1);

fprintf('second delta-v                    %10.4f meters/second \n\n', 1000.0 * dv2);

fprintf('total delta-v                     %10.4f meters/second \n\n\n', 1000.0 * (dv1 + dv2));

fprintf('transfer orbit semimajor axis     %10.4f kilometers \n\n', smat);

fprintf('transfer orbit eccentricity       %10.8f \n\n', ecct);

fprintf('transfer orbit perigee velocity   %10.4f meters/second \n\n', 1000.0 * vt1);

fprintf('transfer orbit apogee velocity    %10.4f meters/second \n\n', 1000.0 * vt2);

fprintf('transfer orbit coast time         %10.4f seconds \n', tof);

fprintf('                                  %10.4f minutes \n', tof / 60.0);

fprintf('                                  %10.4f hours \n\n', tof / 3600.0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create trajectory graphics %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load orbital elements arrays, create state vectors and plot orbits

% Initial orbit
oevi(1) = r1;
oevi(2) = 0.0;
oevi(3) = 0.0;
oevi(4) = 0.0;
oevi(5) = 0.0;
oevi(6) = 0.0;
[ri, vi] = orb2eci(mu, oevi);

%Transfer orbit_first burn
oevti(1) = smat;
oevti(2) = ecct;
oevti(3) = 0.0;
oevti(4) = 0.0;
oevti(5) = 0.0;
oevti(6) = 0.0;
[rti, vti] = orb2eci(mu, oevti);

%Transfer orbit_second burn
oevtf(1) = smat;
oevtf(2) = ecct;
oevtf(3) = 0.0;
oevtf(4) = 0.0;
oevtf(5) = 0.0;
oevtf(6) = 180.0;
[rtf, vtf] = orb2eci(mu, oevtf);

%Final orbit
oevf(1) = r2;
oevf(2) = 0.0;
oevf(3) = 0.0;
oevf(4) = 0.0;
oevf(5) = 0.0;
oevf(6) = 180.0;
[rf, vf] = orb2eci(mu, oevf);

% compute orbital periods

period1 = 2.0 * pi * oevi(1) * sqrt(oevi(1) / mu);  % Period of initial orbit

period2 = 2.0 * pi * oevti(1) * sqrt(oevti(1) / mu);  % Period of transfer orbit

period3 = 2.0 * pi * oevf(1) * sqrt(oevf(1) / mu);   % Period of final orbit

% Divide period for simulation
deltat1 = period1 / 300;
simtime1 = -deltat1;
deltat2 = 0.5 * period2 / 300;
simtime2 = -deltat2;
deltat3 = period3 / 300;
simtime3 = -deltat3;

for i = 1:1:301
    simtime1 = simtime1 + deltat1;
    simtime2 = simtime2 + deltat2;
    simtime3 = simtime3 + deltat3;
    
    % compute initial orbit "normalized" position vector
    [rwrk, vwrk] = twobody2 (mu, simtime1, ri, vi);
    rp1_x(i) = rwrk(1) / req;
    rp1_y(i) = rwrk(2) / req;
    rp1_z(i) = rwrk(3) / req;
    
    % compute transfer orbit position vector
    [rwrk, vwrk] = twobody2 (mu, simtime2, rti, vti);
    rp2_x(i) = rwrk(1) / req;
    rp2_y(i) = rwrk(2) / req;
    rp2_z(i) = rwrk(3) / req;
    
    % compute final orbit posi;tion vector
    [rwrk, vwrk] = twobody2 (mu, simtime3, rf, vf);
    rp3_x(i) = rwrk(1) / req;
    rp3_y(i) = rwrk(2) / req;
    rp3_z(i) = rwrk(3) / req;
end

figure(1);
% create axes vectors
xaxisx = [1 1.5];
xaxisy = [0 0];
xaxisz = [0 0];

yaxisx = [0 0];
yaxisy = [1 1.5];
yaxisz = [0 0];

zaxisx = [0 0];
zaxisy = [0 0];
zaxisz = [1 1.5];

figure (1);
hold on;
grid on;
% plot earth
[x y z] = sphere(24);
h = surf(x, y, z);
colormap([127/255 1 222/255]);
set (h, 'edgecolor', [1 1 1]);

% plot coordinate system axes
plot3(xaxisx, xaxisy, xaxisz, '-g', 'LineWidth', 1);
plot3(yaxisx, yaxisy, yaxisz, '-r', 'LineWidth', 1);
plot3(zaxisx, zaxisy, zaxisz, '-b', 'LineWidth', 1);
% plot initial orbit
plot3(rp1_x, rp1_y, rp1_z, '-r', 'LineWidth', 1.5);
plot3(rp1_x(1), rp1_y(1), rp1_z(1), 'ob');
% plot transfer orbit
plot3(rp2_x, rp2_y, rp2_z, '-b', 'LineWidth', 1.5);
plot3(rp2_x(end), rp2_y(end), rp2_z(end), 'ob');
% plot final orbit
plot3(rp3_x, rp3_y, rp3_z, '-g', 'LineWidth', 1.5);
xlabel('X coordinate (ER)', 'FontSize', 12);
ylabel('Y coordinate (ER)', 'FontSize', 12);
zlabel('Z coordinate (ER)', 'FontSize', 12);
title('Hohmann Transfer: Initial, Transfer and Final Orbits', 'FontSize', 16);

axis equal;
view(50, 20);
rotate3d on;
print -depsc -tiff -r300 hohmann1.eps



