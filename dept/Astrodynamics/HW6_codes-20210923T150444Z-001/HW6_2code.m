r = [6553,2674,3440];
v = [4.326,-1.925,-5.726];

[t,X] = rk4c('cowellformulation',0,86400,1,reshape([r,v],6,1));

r_abs=sqrt(X(1,end)^2+X(2,end)^2+X(3,end)^2)
v_abs=sqrt(X(4,end)^2+X(5,end)^2+X(6,end)^2)

r_answer = X(1:3,end)
v_answer = X(4:6,end)

figure(1)
plot3(X(1,:),X(2,:),X(3,:),'r'), grid on, hold on

figure(2)
subplot(3,2,1)
plot(t,X(1,:),'r'), grid on, hold on, legend('X Position')
subplot(3,2,3)
plot(t,X(2,:),'r'), grid on, hold on, legend('Y Position')
subplot(3,2,5)
plot(t,X(3,:),'r'), grid on, hold on, legend('Z Position')
subplot(3,2,2)
plot(t,X(4,:),'r'), grid on, hold on, legend('X Velocity')
subplot(3,2,4)
plot(t,X(5,:),'r'), grid on, hold on, legend('Y Velocity')
subplot(3,2,6)
plot(t,X(6,:),'r'), grid on, hold on, legend('Z Velocity')