clc
clear all

% 이동할 값 벡터 설정
n = [0 ];

% 주어진 행렬 ff_0
ff_0 = [142 152;
        200 211;
        260 266;
        428 434;
        483 494;
        543 552];

% 그릴 수직선의 x좌표 범위
x_range = 0:719;

% 행의 개수
num_rows = size(ff_0, 1);

% 반복 구간
repeat_range = 720;

% 모든 선이 그려지는 x좌표를 저장할 배열 초기화
all_x_coords = [];

% 플롯 초기화
figure;
hold on;

% 색상 벡터 설정
colors = [1, 0, 0];

% 각 행에 대해 선 그리기
for j = 1:length(n)
    for i = 1:num_rows
        start_x = mod(ff_0(i, 1) + n(j), repeat_range);
        start_x_bf = ff_0(i, 1) + n(j);
        end_x = mod(ff_0(i, 2) + n(j), repeat_range);
        end_x_bf = ff_0(i, 2) + n(j);
       
        % 선 그리기
        if end_x_bf >= 720 && start_x_bf < 720
            plot([start_x ,719], [0, 0], 'Color', [1, 0, 0], 'LineWidth', 2);
            plot([0, end_x_bf - 720], [0, 0], 'Color', [1, 0, 0], 'LineWidth', 2);
            % 해당 구간의 모든 x좌표 저장
            all_x_coords = [all_x_coords, start_x:719, 0:end_x_bf-720];
        elseif start_x_bf >= 720 && end_x_bf >= 720
            plot([start_x_bf - repeat_range, end_x_bf - repeat_range], [0, 0], 'Color', colors(j,:), 'LineWidth', 2);
            % 해당 구간의 모든 x좌표 저장
            all_x_coords = [all_x_coords, start_x_bf - repeat_range:end_x_bf - repeat_range];
        else
            plot([start_x ,end_x], [0, 0], 'Color', [1, 0, 0] , 'LineWidth', 2);
            % 해당 구간의 모든 x좌표 저장
            all_x_coords = [all_x_coords, start_x:end_x];
        end
    end
end

% 플롯 설정
title('영역 연속된 선으로 연결하기');
xlabel('x 좌표');
ylabel('y 좌표');

% 그리드 추가
grid on;

% x축 범위 설정
xlim([0, 719]);

% 플롯 보이기
hold off;

% 모든 선이 그려지는 x좌표 배열을 작업공간에 저장
all_x_coords;

% 중복된 x좌표를 제거하여 유일한 값만 남김
all_unique_x_coords = unique(all_x_coords);

% 작업공간에 유일한 x좌표 저장
all_unique_x_coords;

A = length(all_unique_x_coords) / 720 * 100


length(all_x_coords) - length(all_unique_x_coords) 
