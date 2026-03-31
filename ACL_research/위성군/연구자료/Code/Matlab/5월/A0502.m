clear all

trials = 0;
start_time = tic;
L = 720;
shift = 0;

% while true

% 주어진 행렬 ff_0
ff_00 = ...
[...
22330.89844	22943.67188
32071.34766	34717.44141
43520.97656	45985.66406
55798.71094	57975.11719
67230.23438	69932.46094
78538.06641	80629.27734
...		
];

gg_00 = ...
[...
27972.89063	31660.72266
43305.41016	47868.10547
59775.29297	64328.55469
75999.84375	79516.93359
...		
];


% 저궤도

% spot1 korea

% 22330.89844	22943.67188
% 32071.34766	34717.44141
% 43520.97656	45985.66406
% 55798.71094	57975.11719
% 67230.23438	69932.46094
% 78538.06641	80629.27734

% spot 2 example

% 0	2520
% 11738.78906	13846.28906
% 23932.5	26220.46875
% 35193.63281	37881.79688
% 46760.50781	48004.92188
% 75310.37109	77168.32031
% 85839.02344	0


% 고궤도

% spot 1 korea

% 27972.89063	31660.72266
% 43305.41016	47868.10547
% 59775.29297	64328.55469
% 75999.84375	79516.93359

% spot 2 example

% 0	3810.761719
% 15650.09766	20222.05078
% 32050.48828	36440.44922
% 85431.50391	0


% ff_0 -> n 으로 변환
ff_0 = floor(ff_00/(86400/L));
neg1 = ff_0;

gg_0 = floor(gg_00/(86400/L));
neg2 = gg_0;

% length_ap = ff_0(1,2)-ff_0(1,1)+ff_0(2,2)-ff_0(2,1)+ff_0(3,2)-ff_0(3,1)+...
% +ff_0(4,2)-ff_0(4,1)+4;


% 구간의 시작점과 끝점을 변수에 저장
start_points1 = ff_0(:, 1);
end_points1 = ff_0(:, 2);

start_points2 = gg_0(:, 1);
end_points2 = gg_0(:, 2);

% 모든 구간을 정수로 표시해서 작업 공간에 저장
result = [];
result_1 = [];
result_2 = [];
result_3 = [];
result_4 = [];
result_5 = [];
result_6 = [];
result_7 = [];
result_8 = [];
result_9 = [];
result_10 = [];
result_11 = [];




%% Seed Satellite Setting (0 / 507~512)

L1 = 0;

for n = 0:719
integer_intervals_1 = [];
integer_intervals_2 = [];
integer_intervals_3 = [];

   % 저궤도 seed satellite
   for i = 1:length(ff_0)
       interval_1 = start_points1(i) :end_points1(i) ;
       interval_1(interval_1 > 719) = interval_1(interval_1 > 719) - 720; 
       integer_intervals_1 = [integer_intervals_1, interval_1];
   end

   % 고궤도 seed satellite 틀만 잡기
   for i = 1:length(gg_0)
       interval_2 = start_points2(i) :end_points2(i) ;
       interval_2(interval_2 > 719) = interval_2(interval_2 > 719) - 720; 
       integer_intervals_2 = [integer_intervals_2, interval_2];

   % 고궤도 seed satellite 위치 선정    
       interval_3 = interval_2 + n; 
       interval_3(interval_3 > 719) = interval_3(interval_3 > 719) - 720;
       integer_intervals_3 = [integer_intervals_3, interval_3];
   end

   sum_intervals = ...
       [integer_intervals_1, ...
       integer_intervals_3, ...
       ];

   unique_intervals = unique(sort(sum_intervals));

   difference = length(sum_intervals) - length(unique_intervals) ;

   result_1 = [result_1; n, difference];

end

rows_with_min_1 = result_1(result_1(:, 2) == min(result_1(:, 2)), :);

k_1 = rows_with_min_1(:,1);

random_index_1 = randi([1, length(k_1)]);

random_number_1 = k_1(random_index_1);

A_1 = length(unique_intervals) / 720 * 100;


% renewal

for n = 0:719
integer_intervals_1 = [];
integer_intervals_2 = [];
integer_intervals_3 = [];

   % 저궤도 seed satellite
   for i = 1:length(ff_0)
       interval_1 = start_points1(i) :end_points1(i) ;
       interval_1(interval_1 > 719) = interval_1(interval_1 > 719) - 720; 
       integer_intervals_1 = [integer_intervals_1, interval_1];
   end

   % 고궤도 seed satellite 틀만 잡기
   for i = 1:length(gg_0)
       interval_2 = start_points2(i) :end_points2(i) ;
       interval_2(interval_2 > 719) = interval_2(interval_2 > 719) - 720; 
       integer_intervals_2 = [integer_intervals_2, interval_2];

   % 고궤도 seed satellite 위치 선정    
       interval_3 = interval_2 + random_number_1; 
       interval_3(interval_3 > 719) = interval_3(interval_3 > 719) - 720;
       integer_intervals_3 = [integer_intervals_3, interval_3];
   end

   sum_intervals = ...
       [integer_intervals_1, ...
       integer_intervals_3, ...
       ];

   unique_intervals = unique(sort(sum_intervals));

   difference = length(sum_intervals) - length(unique_intervals) ;

   result_1 = [result_1; n, difference];

end

A_1 = length(unique_intervals) / 720 * 100;



%% 고궤도 위성부터 추가 (L1 / H2)

L1 = 0;

for n = 0:719
integer_intervals_1 = [];
integer_intervals_2 = [];
integer_intervals_3 = [];
integer_intervals_4 = [];

   % 저궤도 seed satellite
   for i = 1:length(ff_0)
       interval_1 = start_points1(i) :end_points1(i) ;
       interval_1(interval_1 > 719) = interval_1(interval_1 > 719) - 720; 
       integer_intervals_1 = [integer_intervals_1, interval_1];
   end

   % 고궤도 seed satellite 틀만 잡기
   for i = 1:length(gg_0)
       interval_2 = start_points2(i) :end_points2(i) ;
       interval_2(interval_2 > 719) = interval_2(interval_2 > 719) - 720; 
       integer_intervals_2 = [integer_intervals_2, interval_2];

   % 고궤도 seed satellite 위치 선정    
       interval_3 = interval_2 + random_number_1; 
       interval_3(interval_3 > 719) = interval_3(interval_3 > 719) - 720;
       integer_intervals_3 = [integer_intervals_3, interval_3];

   % 고궤도 추가   
       interval_4 = interval_2 + n; 
       interval_4(interval_4 > 719) = interval_4(interval_4 > 719) - 720;
       integer_intervals_4 = [integer_intervals_4, interval_4];
   end

   sum_intervals = ...
       [integer_intervals_1, ...
       integer_intervals_3, ...
       integer_intervals_4, ...
       ];

   unique_intervals = unique(sort(sum_intervals));

   difference = length(sum_intervals) - length(unique_intervals) ;

   result_2 = [result_2; n, difference];

end

rows_with_min_2 = result_2(result_2(:, 2) == min(result_2(:, 2)), :);

k_2 = rows_with_min_2(:,1);

random_index_2 = randi([1, length(k_2)]);

random_number_2 = k_2(random_index_2);

A_2 = length(unique_intervals) / 720 * 100;


% renewal

for n = 0:719
integer_intervals_1 = [];
integer_intervals_2 = [];
integer_intervals_3 = [];
integer_intervals_4 = [];

   % 저궤도 seed satellite
   for i = 1:length(ff_0)
       interval_1 = start_points1(i) :end_points1(i) ;
       interval_1(interval_1 > 719) = interval_1(interval_1 > 719) - 720; 
       integer_intervals_1 = [integer_intervals_1, interval_1];
   end

   % 고궤도 seed satellite 틀만 잡기
   for i = 1:length(gg_0)
       interval_2 = start_points2(i) :end_points2(i) ;
       interval_2(interval_2 > 719) = interval_2(interval_2 > 719) - 720; 
       integer_intervals_2 = [integer_intervals_2, interval_2];

   % 고궤도 추가   
       interval_4 = interval_2 + random_number_2; 
       interval_4(interval_4 > 719) = interval_4(interval_4 > 719) - 720;
       integer_intervals_4 = [integer_intervals_4, interval_4];
   end



   sum_intervals = ...
       [integer_intervals_1, ...
       integer_intervals_3, ...
       integer_intervals_4, ...
       ];

   unique_intervals = unique(sort(sum_intervals));

   difference = length(sum_intervals) - length(unique_intervals) ;

   result_2 = [result_2; n, difference];

end

A_2 = length(unique_intervals) / 720 * 100;




%% Acess Profile
% x축 범위 설정
x_range = 0:L-1;

figure;
subplot(3, 1, 1);

hex_color1 = '#005eb5';

red1 = hex2dec(hex_color1(2:3)) / 255; % R 값
green1 = hex2dec(hex_color1(4:5)) / 255; % G 값
blue1 = hex2dec(hex_color1(6:7)) / 255; % B 값

hex_color11 = '#6ca6d6';

red11 = hex2dec(hex_color11(2:3)) / 255; % R 값
green11 = hex2dec(hex_color11(4:5)) / 255; % G 값
blue11 = hex2dec(hex_color11(6:7)) / 255; % B 값


% 각 구간별로 y=1인 위치 표시
for i = 1:size(ff_0, 1)
    if neg1(i) < 0
        start_time_new = L + neg1(i,1);
        end_time_new = L + neg1(i,2);
    else 
        start_time_new = ff_0(i, 1);
        end_time_new = ff_0(i, 2);
    end
    start_time = ff_0(i, 1);
    end_time = ff_0(i, 2);

    if end_time_new > L-1
        end_time_new = L-1;
    end

    % 시작 시간과 종료 시간에 해당하는 x 인덱스 찾기
    start_index = find(x_range >= start_time, 1);
    end_index = find(x_range <= end_time, 1, 'last');
    stn = start_time_new;
    etn = end_time_new;

    % y 값 설정
    y_values1a = zeros(size(x_range));
    y_values1a(start_index:end_index) = 1;
    if neg1(i) < 0
        y_values1a(stn:etn) = 1;
    end

    % 그래프 그리기
    plot(x_range, y_values1a, '-o r', 'LineWidth', 1, 'MarkerSize', 1, 'Color', [red1, green1, blue1]);
    hold on;
    area(x_range, y_values1a, 'FaceColor', [red11, green11, blue11]);
    plot(x_range, y_values1a, 'Color', [red1, green1, blue1], 'LineWidth', 1.5); 
end

% 그래프 세부 설정
xlabel('n');
ylabel('v0j[n]');
title('Access Profile');
grid on;
xlim([0, L-1]);
ylim([0, 1]);

%% Constellation Pattern Vector

x_range = 0:L-1;

% 720x2 행렬 초기화
sorted_matrix1 = zeros(720, 2);
sorted_matrix1(:, 1) = 0:719;
sorted_matrix1(integer_intervals_1 + 1, 2) = 1; 

% n의 값
n_L = [random_number_1] - 1;

% 새로운 행렬을 생성하여 값을 이동
shifted_matrix = zeros(720, 2);

% 1열은 기존 행렬의 1열과 동일
shifted_matrix(:, 1) = sorted_matrix1(:, 1);

% 2열의 위치를 n만큼 이동하고, 720을 초과하는 경우 0부터 다시 시작
nonzero_indices = find(sorted_matrix1(:, 2) == 1);  % 1인 위치 찾기
shifted_indices = mod(nonzero_indices + n_L, 720);  % 720을 초과할 경우 순환

% 새 위치에 1 설정
shifted_matrix(shifted_indices + 1, 2) = 1;  % MATLAB의 1 기반 인덱싱 고려


subplot(3, 1, 2);

    % y 값 설정
    y_values1p = zeros(size(x_range));
    y_values1p(random_number_1(1)) = 1;


plot(x_range, y_values1p, 'Color', [red1, green1, blue1], 'LineWidth', 1, 'MarkerSize', 1);
title('Constellation Pattern Vector');
hold on;

xlim([0, L-1]);
ylim([0, 1]);

% 그래프 세부 설정
xlabel('n');
ylabel('x[n]');
grid on;
xlim([0, L-1]);
ylim([0, 1]);


%% Coverage Timline

x_range = 0:L-1;

subplot(3, 1, 3);
plot(shifted_matrix(:,1), shifted_matrix(:,2), 'Color', [red1, green1, blue1], 'LineWidth', 1, 'MarkerSize', 1);
title('Coverage Timeline');
hold on;
area(shifted_matrix(:,1), shifted_matrix(:,2), 'FaceColor', [red11, green11, blue11]);
plot(shifted_matrix(:,1), shifted_matrix(:,2), 'Color', [red1, green1, blue1], 'LineWidth', 1);

xlim([0, L-1]);
ylim([0, 3]);

% 그래프 세부 설정
xlabel('n');
ylabel('bj[n]');
grid on;
xlim([0, L-1]);



% random_N(1) = (random_number_1);
% random_N(2) = (random_number_2);
% random_N(3) = (random_number_3);
% random_N(4) = (random_number_4);
% random_N(5) = (random_number_5);
% random_N(6) = (random_number_6);
% random_N(7) = (random_number_7);
% random_N(8) = (random_number_8)
% 
%   if A_8 >= 99
% 
% 
%         break;
%     end
%  trials = trials + 1; 
% end
% 
% % 시뮬레이션 종료 시간 기록
% elapsed_time = toc(start_time);
% 
% % 결과 출력
% fprintf('시행 횟수: %d\n', trials);
% fprintf('걸린 시간(초): %.2f\n', elapsed_time);

