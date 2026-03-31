% x값 생성 (0부터 720까지 1씩 증가)
x = 0:720;

% 원하는 x 값들
desired_x_values = [0, 33, 65, 98, 131, 164, 196, 229, 262, ...
    295, 327, 360, 393, 425, 458, 491, 524, 556, 589, 622, 655, 687];

% y값 생성 (0과 1 사이의 값)
y = zeros(size(x)); % 일단 모든 값은 0으로 초기화

% 원하는 x 값에 대응하는 y 값을 1로 설정
for i = 1:length(desired_x_values)
    x_val = desired_x_values(i);
    y(x == x_val) = 1;
end

% 그래프 그리기
plot(x, y, 'r');
xlim([0 720]); % x축 범위를 0부터 720까지로 제한
xlabel('n');
ylabel('x[n]');
title('Constellation pattern vector');
