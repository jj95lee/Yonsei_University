% 2019-10-21
%case 1a
A = [ 0     1     0     0;
    -1    -1     0    -1;
    -1     0     0     1;
     0    -1    -1     1];

exA=expm(A);

% exA =
% 
%     0.5858    0.6524    0.1663   -0.4972
%    -0.8187    0.4306    0.4972   -0.9833
%    -0.4845   -0.8282    0.2533    1.8131
%     1.1607   -0.8171   -1.6468    2.3973

temp1=exA(1:2,1:2);

% temp1 =
% 
%     0.5858    0.6524
%    -0.8187    0.4306

temp2=exA(1:2,3:4);

% temp2 =
% 
%     0.1663   -0.4972
%     0.4972   -0.9833

lambda0=-inv(temp2)*temp1*[1;0];

% ans =
% 
%    11.7414
%     5.1045


for i=1:101
result(:,i)=expm(A*(i-1)/100)*[1;0;lambda0];
end
t=0:0.01:1;
figure()
plot(t,result(1,:),'linewidth',1.6);% x1
hold on;
plot(t,result(2,:),'linewidth',1.6);% x2
plot(t,-result(4,:),'linewidth',1.6);% u
legend({'x1(t)','x1(t)','u(t)'})
% plot(simout.time,simout.data(:,1),'linewidth',1.6)
xlabel('time (sec)');
title('case 1a')
% xlim([0:1]);
grid on;
set(gca,'fontsize',16);

% case 1b

x0=[1;1];
 a=exA(2:3,1:2)

% a =
% 
%    -0.8187    0.4306
%    -0.4845   -0.8282

b=exA(2:3,3:4)

% b =
% 
%     0.4972   -0.9833
%     0.2533    1.8131

lambda0b=-inv(b)*a*[1;1];

% ans =
% 
%     1.7334
%     0.4819

for i=1:101
resultb(:,i)=expm(A*(i-1)/100)*[1;1;lambda0b];
end
t=0:0.01:1;
figure()
plot(t,resultb(1,:),'linewidth',1.6);% x1
hold on;
plot(t,resultb(2,:),'linewidth',1.6);% x2
plot(t,-resultb(4,:),'linewidth',1.6);% u
legend({'x1(t)','x1(t)','u(t)'})
xlabel('time (sec)');
title('case 1b')
grid on;
set(gca,'fontsize',16);