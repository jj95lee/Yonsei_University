%% Problem 1(a)
wx = 0.1; % [rad s^-1]
wy = -0.05; % [rad s^-1]
wz = 0.03; % [rad s^-1]


%% Problem 1(a) - (i)
Sw = [0     -wz     wy  ;
      wz    0       -wx ;
      -wy   wx      0  ];
answer_1ai = expm(Sw*10);


%% Problem 1(a) - (ii)
answer_1aii = [0; 0; 0; acos( (trace(answer_1ai) - 1)/2 )];
answer_1aii(1) = (answer_1ai(3, 2) - answer_1ai(2, 3))/2/sin(answer_1aii(4));
answer_1aii(2) = (answer_1ai(1, 3) - answer_1ai(3, 1))/2/sin(answer_1aii(4));
answer_1aii(3) = (answer_1ai(2, 1) - answer_1ai(1, 2))/2/sin(answer_1aii(4));
answer_1aii(1:3) * answer_1aii(4)


%% Problem 1(a) - (iii)
syms phi theta psi
answer_1aiii_sym = [1 0 0; 0 cos(phi) sin(phi); 0 -sin(phi) cos(phi)] * ...
                   [cos(theta) 0 -sin(theta); 0 1 0; sin(theta) 0 cos(theta)] *...
                   [cos(psi) sin(psi) 0; -sin(psi) cos(psi) 0; 0 0 1];
theta = rad2deg(asin(-answer_1ai(1, 3)))
psi = rad2deg(asin(answer_1ai(1, 2)/cosd(theta)))
phi = rad2deg(asin(answer_1ai(2, 3)/cosd(theta)))
               

%% Problem 1(a) - (iv)
q = [answer_1aii(1:3) * sin(answer_1aii(4)/2)];
q(4, 1) = cos(answer_1aii(4)/2);


%% Problem 1(b) - (ii)
answer_1bii = answer_1ai*[wx; wy; wz];


%% Video
f311 = figure(311);
set(f311, 'Position', [100 100 700 700])

for idx = 1:101
    
    % t = (idx - 1)/10;
    rm = expm(Sw*(idx - 1)/10); % Rotation matrix(t)
    
    hold on
    axis square
    
    quiver3(-1, 0, 0, 2, 0, 0, 'Color', 'black', 'LineWidth', 1)
    quiver3(0, -1, 0, 0, 2, 0, 'Color', 'black', 'LineWidth', 1)
    quiver3(0, 0, -1, 0, 0, 2, 'Color', 'black', 'LineWidth', 1)
    
    text(1, 0, 0, 'I', 'FontWeight', 'bold', 'FontSize', 15)
    text(0, 1, 0, 'J', 'FontWeight', 'bold', 'FontSize', 15)
    text(0, 0, 1, 'K', 'FontWeight', 'bold', 'FontSize', 15)
    
    axis([-1 1 -1 1 -1 1])
    xlabel('I'); ylabel('J'); zlabel('K')
    
    % Angular velocity
    w = rm*[0.1; -0.05; 0.03]; % Actual angular velocity
    W = w/norm(w); % Angular velocity normalized
    quiver3(0, 0, 0, W(1), W(2), W(3), 'Color', 'c', 'LineWidth', 2) % Angular velocity vector
    text(W(1), W(2), W(3), ['I\omegaI = ', num2str(norm(w))], 'FontWeight', 'bold', 'FontSize', 10) % Name
    
    % Angle
    x = rm*[1; 0; 0]; % i-axis in inertial frame
    y = rm*[0; 1; 0]; % j-axis in inertial frame
    z = rm*[0; 0; 1]; % k-axis in inertial frame
    
    quiver3(0, 0, 0, x(1), x(2), x(3), 'Color', [0.8510 0.3294 0.1020], 'LineWidth', 3)
    quiver3(0, 0, 0, y(1), y(2), y(3), 'Color', [0.4706 0.6706 0.1882], 'LineWidth', 3)
    quiver3(0, 0, 0, z(1), z(2), z(3), 'Color', [0.0000 0.4510 0.7412], 'LineWidth', 3)
    
    text(x(1), x(2), x(3), 'i', 'FontWeight', 'bold', 'FontSize', 10, 'FontName', 'Calibri')
    text(y(1), y(2), y(3), 'j', 'FontWeight', 'bold', 'FontSize', 10, 'FontName', 'Calibri')
    text(z(1), z(2), z(3), 'k', 'FontWeight', 'bold', 'FontSize', 10, 'FontName', 'Calibri')
    
    view(135, 30)
    grid on
    title(['t = ' num2str((idx - 1)/10) ' sec'])
    
    filename = ['HW3_1_frame' num2str(idx, '%03.f') '.png'];
    saveas(gcf, filename);
    clf
end
close all