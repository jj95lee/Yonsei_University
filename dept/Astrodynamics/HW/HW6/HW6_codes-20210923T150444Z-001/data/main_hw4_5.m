xi = [6553;2674;3440;4.326;-1.925;-5.726];  % for HW4

opt = 2;    % 1:two-body / 2:j2,j3,j4 / 3:all perturbation

options = odeset('RelTol',1e-12,'AbsTol',1e-12);
[T Y] = ode45(@(t,x) hw4_5(t,x,opt),[0:10:0.2*86400],xi,options);

figure(1)
plot3(Y(:,1),Y(:,2),Y(:,3),'r'), grid on, hold on

figure(2)
subplot(3,2,1)
plot(T/86400,Y(:,1),'r'), grid on, hold on, legend('X Position')
subplot(3,2,3)
plot(T/86400,Y(:,2),'r'), grid on, hold on, legend('Y Position')
subplot(3,2,5)
plot(T/86400,Y(:,3),'r'), grid on, hold on, legend('Z Position')
subplot(3,2,2)
plot(T/86400,Y(:,4),'r'), grid on, hold on, legend('X Velocity')
subplot(3,2,4)
plot(T/86400,Y(:,5),'r'), grid on, hold on, legend('Y Velocity')
subplot(3,2,6)
plot(T/86400,Y(:,6),'r'), grid on, hold on, legend('Z Velocity')