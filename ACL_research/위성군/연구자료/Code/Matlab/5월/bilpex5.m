clear all

%% n의 배열
n =  [65, 144, 285, 361];



NumberofSatellite = length(n);

%% freefleyer를 통해 구한 seed satellite의 pass 구간

ff_0 = ...
[12 32
106	126
199	216
293	299
462	468
544	562
634	655]-360;

% 시작 시각과 종료 시각이 저장된 행렬
start_time = ff_0(:, 1);
end_time = ff_0(:, 2);


%ff_0 = seedsatellite
ff_1 = mod(ff_0 + n(1),720);
ff_2 = mod(ff_0 + n(2),720);
ff_3 = mod(ff_0 + n(3),720);
ff_4 = mod(ff_0 + n(4),720);
% ff_5 = ff_0 + n(5);
% ff_6 = ff_0 + n(6);
% ff_7 = ff_0 + n(7);
% ff_8 = ff_0 + n(8);
% ff_9 = ff_0 + n(9);
% ff_10 = ff_0 + n(10);
% ff_11 = ff_0 + n(11);
% ff_12 = ff_0 + n(12);
% ff_13 = ff_0 + n(13);
% ff_14 = ff_0 + n(14);
% ff_15 = ff_0 + n(15);
% ff_16 = ff_0 + n(16);
% ff_17 = ff_0 + n(17);
% ff_18 = ff_0 + n(18);
% ff_19 = ff_0 + n(19);
% ff_20 = ff_0 + n(20);
% ff_21 = ff_0 + n(21);
% ff_22 = ff_0 + n(22);
% ff_23 = ff_0 + n(23);
% ff_24 = ff_0 + n(24);


% x축의 범위 설정 (0부터 716까지)
x_range =  0:1:716;
%% 여기부터
% 그래프 초기화
figure;   


% x축의 범위 설정 (0부터 86400까지)
x_range = 0:716;


%% 모든 작업공간 그래프 그리기
for idx = 1:NumberofSatellite
    ff = eval(['ff_', num2str(idx)]);

    for i = 1:size(ff,1)
        start_time = ff(i, 1);
        end_time = ff(i, 2);

        start_index = find(x_range >= start_time, 1);
        end_index = find(x_range <= end_time, 1, 'last');

        y_range = zeros(size(x_range));
        y_range(start_index:end_index) = 1;


% 첫 번째 구간 추가
range1 = x_range(start_index:end_index);

% 두 번째 구간 추가 (만약 end_time이 717을 초과하는 경우)
if end_time > 716 & start_time < 717
    range2 = [x_range(1:end_time-716)];
elseif start_time > 716
    range2 = [x_range(start_time-716:end_time-717)];
 else
    range2 = range1;
end


% 그래프 그리기
 fill([x_range(start_index:end_index), fliplr(x_range(start_index:end_index))], ...
             [y_range(start_index:end_index), zeros(1, length(start_index:end_index))], ...
             [0.3, 0.6, 0.9], 'EdgeColor', 'none');
        hold on;

 fill([range2, fliplr(range2)], ...
             [ones(1, length(range2)), zeros(1, length(range2))], ...
             [0.3, 0.6, 0.9]);

 xlim([0, 716]);
    end
end

xlabel("n");  % 초 단위로 x축 표시
ylabel("v0j[n]");
title('Coverage Timeline')

%%
%% 여기전까지 
for idx = 1:NumberofSatellite
    ff = eval(['ff_', num2str(idx)]);
    x_coords{idx} = []; % 해당 작업공간의 x좌표를 저장할 배열 초기화

    for i = 1:size(ff,1)
        start_time = ff(i, 1);
        end_time = ff(i, 2);

        % y=1일 때의 시작 시간과 종료 시간의 x좌표 계산
        if end_time > 716 & start_time < 717
           start_index = find(x_range >= start_time, 1);
           end_index = find(x_range < end_time-717, 1, 'last');
        elseif start_time > 716
           start_index = find(x_range >= start_time-717, 1);
           end_index = find(x_range <= end_time-717, 1, 'last');
        else
        start_index = find(x_range >= start_time, 1);
        end_index = find(x_range <= end_time, 1, 'last');
        end



        % 시작 시간부터 종료 시간까지의 x좌표를 배열에 저장
        if end_time > 716 & start_time < 717
          x_coords{idx} = [x_coords{idx}, x_range(start_index:717),(0:end_index)];

        else
          x_coords{idx} = [x_coords{idx}, x_range(start_index:end_index)];
        end



    end               
end
%
% 모든 작업공간의 x 좌표를 모은 배열 초기화
all_x_coords = [];

%% 각 작업공간의 x 좌표를 all_x_coords에 추가하고 작은 수부터 오름차순으로 정렬
for idx = 1:NumberofSatellite
    all_x_coords = [all_x_coords, x_coords{idx}];
end

all_x_coords = sort(all_x_coords);

% 중복된 값을 제거하고 유일한 값들만 추출
unique_x_coords = unique(all_x_coords);

coverage1 = length(unique_x_coords) / 717 * 100;
% 
% 
% % 그래프 세부 설정
% xlabel('시간 (초)');
% ylabel('Y값');
% title('구간 그래프');
% grid on;
% xlim([0, 716]);
% 
% 
% %%
% 
% % 초기화
% RV = [];
% current_value = all_x_coords(1);
% start_index = 1;
% 
% % 각 원소에 대해 반복하면서 연속된 값의 구간을 찾음
% for i = 2:length(all_x_coords)
%     if all_x_coords(i) ~= current_value
%         % 현재 값과 다른 값이 나타난 경우, 해당 구간의 시작과 끝을 저장하고 현재 값 갱신
%         end_index = i - 1;
%         RV = [RV; current_value, end_index - start_index + 1];
%         current_value = all_x_coords(i);
%         start_index = i;
%     end
% end
% 
% % 마지막 구간에 대한 처리
% end_index = length(all_x_coords);
% RV = [RV; current_value, end_index - start_index + 1];
% 
% 
% %% 빠진 숫자 찾기
% 
% % 주어진 행렬
% matrix = RV ;
% 
% % 현재 1열에서 빠진 숫자를 찾기
% missing_numbers = setdiff(0:716, matrix(:,1));
% 
% % 새롭게 추가된 행에 대해 0으로 입력
% new_rows = zeros(length(missing_numbers), 2);
% new_rows(:, 1) = missing_numbers;
% 
% % 새로운 행렬 생성
% new_matrix = [matrix; new_rows];
% 
% % 1열을 기준으로 행렬을 정렬
% sorted_matrix = sortrows(new_matrix, 1);
% 
% 
% %% Acess Profile
% % x축 범위 설정
% x_range = 0:716;
% 
% figure;
% subplot(3, 1, 1);
% 
% % 각 구간별로 y=1인 위치 표시
% for i = 1:size(ff_0, 1)
%     start_time = ff_0(i, 1);
%     end_time = ff_0(i, 2);
% 
%     % 시작 시간과 종료 시간에 해당하는 x 인덱스 찾기
%     start_index = find(x_range >= start_time, 1);
%     end_index = find(x_range <= end_time, 1, 'last');
% 
%     % y 값 설정
%     y_values = zeros(size(x_range));
%     y_values(start_index:end_index) = 1;
% 
%     % 그래프 그리기
%     plot(x_range, y_values, '-o r', 'LineWidth', 1, 'MarkerSize', 1, 'Color', [0, 0, 0.3]);
%     hold on;
% end
% 
% % 그래프 세부 설정
% xlabel('n');
% ylabel('v0j[n]');
% title('Access Profile');
% grid on;
% xlim([0, 716]);
% ylim([0, 1]);
% 
% %% Constellation Pattern Vector
% 
% x_range = 0:716;
% 
% subplot(3, 1, 2);
% 
%     % y 값 설정
%     y_values = zeros(size(x_range));
%     y_values(n(1,:)) = 1;
% 
% 
% plot(x_range, y_values, 'Color', [1, 0, 0], 'LineWidth', 1, 'MarkerSize', 1);
% title('Constellation Pattern Vector');
% xlim([0, 716]);
% ylim([0, 1]);
% 
% % 그래프 세부 설정
% xlabel('n');
% ylabel('x[n]');
% grid on;
% xlim([0, 716]);
% ylim([0, 1]);
% 
% 
% %% Coverage Timline
% 
% x_range = 0:716;
% 
% subplot(3, 1, 3);
% plot(sorted_matrix(:,1), sorted_matrix(:,2), 'Color', [0, 0.5, 0], 'LineWidth', 1, 'MarkerSize', 1);
% title('Coverage Timeline');
% xlim([0, 716]);
% ylim([0, 4]);
% 
% % 그래프 세부 설정
% xlabel('n');
% ylabel('bj[n]');
% grid on;
% xlim([0, 716]);


% %% V
% % 작업공간을 저장할 변수 초기화
% acp = zeros(1, 717);
% 
% % 각 구간에 대해 1로 설정
% for i = 1:size(ff_0, 1)
%     start_index = ff_0(i, 1);
%     end_index = ff_0(i, 2);
%     acp(start_index:end_index) = 1;
% end
% 
% 
% shiftAmount = 1;
% 
% V = [];
% 
% for i = 1:length(acp)
%     newRow = circshift(acp, [0, i-shiftAmount]);
%     V = [V; newRow];
% end
% 
% %% Z
% 
% % 행렬 크기 정의
% row_count = 717;
% col_count = 1;
% 
% % 모두 0으로 채워진 행렬 생성
% Z = zeros(row_count, col_count);
% 
% % 덮어쓸 특정 행 지정
% % target_rows = [39, 73, 79, 89, 170, 184, 234, 250, 331, 341, 347, 492, 502, 542, 638, 648, 654, 663];
% % target_rows = [18, 114, 124, 130, 139, 235, 269, 275, 285, 366, 380, 430, 446, 527, 537, 543, 688, 698];
% 
% % 특정 행에 1 덮어쓰기
% Z(n, :) = 1;
% 
% %% X Y
% 
% X = V*Z;
% 
% Y = sum(X >= 1)/717*100;
% 
% % 원하는 1의 개수 설정
% desired_ones = 18;
% 
% % 임계값 계산
% threshold = desired_ones / 717;
% 
% % 시행 횟수 초기화
% trials = 0;
% 
% % 시뮬레이션 시작 시간 기록
% start_time = tic;
% 
% % 결과가 나올 때까지 반복
% while true
%     % B 행렬을 랜덤하게 생성
%     B = rand(717, 1) < threshold;
% 
%     % 생성된 B 행렬에서 1의 개수가 desired_ones 보다 클 때까지 반복
%     while sum(B) > desired_ones
%         B = rand(717, 1) < threshold;
%         trials = trials + 1; % 시행 횟수 증가
%     end
% 
%     % V와 B 행렬의 곱셈을 통해 C 행렬을 계산
%     C = V * B;
% 
%     % C 행렬의 조건 확인
%     if all(sum(C(:) >= 1) / 717 >= 0.93)
%         break; % 조건을 만족하면 반복 중단
%     end
% 
%     trials = trials + 1; % 시행 횟수 증가
% end
% 
% % 시뮬레이션 종료 시간 기록
% elapsed_time = toc(start_time);
% 
% % 결과 출력
% fprintf('시행 횟수: %d\n', trials);
% fprintf('걸린 시간(초): %.2f\n', elapsed_time);
% 
% cov_percent = sum(C(:) >= 1) / 717 * 100

duplicates = sum(diff(sort(all_x_coords)) == 0);
