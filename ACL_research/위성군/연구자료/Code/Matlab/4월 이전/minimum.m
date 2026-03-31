
% 주어진 A 행렬
% A = [1, 1, 1, 0, 0, 1, 1, 0, 0, 1;
%      1, 1, 1, 1, 0, 0, 1, 1, 0, 0;
%      0, 1, 1, 1, 1, 0, 0, 1, 1, 0;
%      0, 0, 1, 1, 1, 1, 0, 0, 1, 1;
%      1, 0, 0, 1, 1, 1, 1, 0, 0, 1;
%      1, 1, 0, 0, 1, 1, 1, 1, 0, 0;
%      0, 1, 1, 0, 0, 1, 1, 1, 1, 0;
%      0, 0, 1, 1, 0, 0, 1, 1, 1, 1;
%      1, 0, 0, 1, 1, 0, 0, 1, 1, 1;
%      1, 1, 0, 0, 1, 1, 0, 0, 1, 1];

% 원하는 1의 개수 설정
desired_ones = 34;

% 임계값 계산
threshold = desired_ones / 720;

% C 행렬의 조건과 1의 개수를 충족하는 B 행렬 찾기
C_condition = false;

while ~C_condition
    % 1의 확률을 높게 설정
    B = rand(720, 1) < threshold;
    while sum(B) > desired_ones % 1의 개수가 N 이하일 때까지 반복
        B = rand(720, 1) < threshold;
    end
    C = V * B; % C 행렬 계산
    C_condition = all(C >= 1); % C 행렬의 조건 확인
end

% B 행렬 그래프
subplot(2, 1, 1);
plot(B, '-o r', 'LineWidth', 1, 'MarkerSize', 1);
title('Constellation Pattern Vector');
xlim([0, 720]);
ylim([-0.2, 1.5]);

% C 행렬 그래프
subplot(2, 1, 2);
plot(find(C), C(C~=0), '-o', 'Color', [0, 0.5, 0], 'LineWidth', 1, 'MarkerSize', 1);
title('Coverage Timeline');
xlim([0, 720]);
ylim([0, 10]);

xlabel('n');
ylabel('');