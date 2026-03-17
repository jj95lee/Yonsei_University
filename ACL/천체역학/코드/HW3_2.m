%% Problem 2(b)
temp_stateeq = @(t, state) HW3_stateeq(t, state, [2000, 2000, 1000]);
[list_t2b, list_state2b] = ode45(temp_stateeq, 0:0.1:1000, [0, 0, 0.8, 0, 0, 0]);

f1 = figure(321);
set(f1, 'Position', [100 100 1200 800])
subplot(3, 2, 2)
plot(list_t2b, list_state2b(:, 4));
xlabel('t (sec)'); ylabel('\psi (rad)');
grid on

subplot(3, 2, 4)
plot(list_t2b, list_state2b(:, 5));
xlabel('t (sec)'); ylabel('\theta (rad)');
grid on

subplot(3, 2, 6)
plot(list_t2b, list_state2b(:, 6));
xlabel('t (sec)'); ylabel('\phi (rad)');
grid on

subplot(3, 2, 1)
plot(list_t2b, list_state2b(:, 1));
xlabel('t (sec)'); ylabel('\omega_x (rad sec^{-1})');
grid on

subplot(3, 2, 3)
plot(list_t2b, list_state2b(:, 2));
xlabel('t (sec)'); ylabel('\omega_y (rad sec^{-1})');
grid on

subplot(3, 2, 5)
plot(list_t2b, list_state2b(:, 3));
xlabel('t (sec)'); ylabel('\omega_z (rad sec^{-1})');
grid on

close all


%% Video
psi = list_state2b(:, 4);
theta = list_state2b(:, 5);
phi = list_state2b(:, 6);

f321 = figure(321);
set(f321, 'Position', [100 100 700 700])

for idx = 1:101
    
    % t = (idx - 1)/10;
    % Form rotation matrix
    Rx = [1 0 0; 0 cos(phi(idx)) sin(phi(idx)); 0 -sin(phi(idx)) cos(phi(idx))];
    Ry = [cos(theta(idx)) 0 -sin(theta(idx)); 0 1 0; sin(theta(idx)) 0 cos(theta(idx))];
    Rz = [cos(psi(idx)) sin(psi(idx)) 0; -sin(psi(idx)) cos(psi(idx)) 0 ; 0 0 1];
    rm = Rx*Ry*Rz;
    
    hold on
    axis square
    
    quiver3(-1, 0, 0, 2, 0, 0, 'Color', 'black', 'LineWidth', 2)
    quiver3(0, -1, 0, 0, 2, 0, 'Color', 'black', 'LineWidth', 2)
    quiver3(0, 0, -1, 0, 0, 2, 'Color', 'black', 'LineWidth', 2)
    
    text(1, 0, 0, 'I', 'FontWeight', 'bold', 'FontSize', 13)
    text(0, 1, 0, 'J', 'FontWeight', 'bold', 'FontSize', 13)
    text(0, 0, 1, 'K', 'FontWeight', 'bold', 'FontSize', 13)
    
    axis([-1 1 -1 1 -1 1])
    xlabel('I'); ylabel('J'); zlabel('K')
    
    % Angular velocity
    w = rm*list_state2b(idx, 1:3)'; % Actual angular velocity
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
    
    filename = ['HW3_2b_frame' num2str(idx, '%03.f') '.png'];
    saveas(gcf, filename);
    clf
end
close all

%% Problem 2(c)
[list_t2c, list_state2c] = ode45(temp_stateeq, 0:0.1:1000, [0.05, -0.05, 0, 0, 0, 0]);

f1 = figure(322);
set(f1, 'Position', [100 100 1200 800])
subplot(3, 2, 2)
plot(list_t2c, list_state2c(:, 4));
xlabel('t (sec)'); ylabel('\psi (rad)');
grid on

subplot(3, 2, 4)
plot(list_t2c, list_state2c(:, 5));
xlabel('t (sec)'); ylabel('\theta (rad)');
grid on

subplot(3, 2, 6)
plot(list_t2c, list_state2c(:, 6));
xlabel('t (sec)'); ylabel('\phi (rad)');
grid on

subplot(3, 2, 1)
plot(list_t2c, list_state2c(:, 1));
xlabel('t (sec)'); ylabel('\omega_x (rad sec^{-1})');
grid on

subplot(3, 2, 3)
plot(list_t2c, list_state2c(:, 2));
xlabel('t (sec)'); ylabel('\omega_y (rad sec^{-1})');
grid on

subplot(3, 2, 5)
plot(list_t2c, list_state2c(:, 3));
xlabel('t (sec)'); ylabel('\omega_z (rad sec^{-1})');
grid on

close all


%% Video
psi = list_state2c(:, 4);
theta = list_state2c(:, 5);
phi = list_state2c(:, 6);

f322 = figure(322);
set(f322, 'Position', [100 100 700 700])

for idx = 1:101
    
    % t = (idx - 1)/10;
    % Form rotation matrix
    Rx = [1 0 0; 0 cos(phi(idx)) sin(phi(idx)); 0 -sin(phi(idx)) cos(phi(idx))];
    Ry = [cos(theta(idx)) 0 -sin(theta(idx)); 0 1 0; sin(theta(idx)) 0 cos(theta(idx))];
    Rz = [cos(psi(idx)) sin(psi(idx)) 0; -sin(psi(idx)) cos(psi(idx)) 0 ; 0 0 1];
    rm = Rx*Ry*Rz;
    
    hold on
    axis square
    
    quiver3(-1, 0, 0, 2, 0, 0, 'Color', 'black', 'LineWidth', 2)
    quiver3(0, -1, 0, 0, 2, 0, 'Color', 'black', 'LineWidth', 2)
    quiver3(0, 0, -1, 0, 0, 2, 'Color', 'black', 'LineWidth', 2)
    
    text(1, 0, 0, 'I', 'FontWeight', 'bold', 'FontSize', 13)
    text(0, 1, 0, 'J', 'FontWeight', 'bold', 'FontSize', 13)
    text(0, 0, 1, 'K', 'FontWeight', 'bold', 'FontSize', 13)
    
    axis([-1 1 -1 1 -1 1])
    xlabel('I'); ylabel('J'); zlabel('K')
    
    % Angular velocity
    w = rm*list_state2c(idx, 1:3)'; % Actual angular velocity
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
    
    filename = ['HW3_2c_frame' num2str(idx, '%03.f') '.png'];
    saveas(gcf, filename);
    clf
end
close all