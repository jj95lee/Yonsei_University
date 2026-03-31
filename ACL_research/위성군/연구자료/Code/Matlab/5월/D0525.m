clc
clear all

trials = 0;
start_time = tic;
L = 100;

while true

% pyeongyang

ff_0 = ...
[...
0 5
25 30
40 45
50 55
80 85
...
];


% length_ap = ff_0(1,2)-ff_0(1,1)+ff_0(2,2)-ff_0(2,1)+ff_0(3,2)-ff_0(3,1)+...
% +ff_0(4,2)-ff_0(4,1)+ff_0(5,2)-ff_0(5,1)+ff_0(6,2)-ff_0(6,1)+ff_0(7,2)-ff_0(7,1)+7;


% 구간의 시작점과 끝점을 변수에 저장
start_points1 = ff_0(:, 1);
end_points1 = ff_0(:, 2);

% 구간의 시작점과 끝점을 변수에 저장
start_points = ff_0(:, 1);
end_points = ff_0(:, 2);


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


%% Random35 (위성 N대)

integer_intervals_1 = [];
integer_intervals_2 = [];
integer_intervals_3 = [];
integer_intervals_4 = [];
integer_intervals_5 = [];
% integer_intervals_6 = [];
% integer_intervals_7 = [];
% integer_intervals_8 = [];
% integer_intervals_9 = [];
% integer_intervals_10 = [];

% 변수 목록
vars = who;

% 'integer_intervals_'로 시작하는 변수들만 필터링
pattern = '^integer_intervals_\d+$'; % 정규 표현식 패턴
N = 0;

for i = 1:length(vars)
    if ~isempty(regexp(vars{i}, pattern, 'once'))
        N = N + 1;
    end
end


  for i = 1:length(ff_0)
      % satellite 1 위치(seed)
      interval_1 = start_points(i) :end_points(i) ;
      interval_1(interval_1 > (L-1)) = interval_1(interval_1 > (L-1)) - L; 
      integer_intervals_1 = [integer_intervals_1, interval_1];

      % satellite 2 위치
      interval_2 = interval_1 + round(L/N); 
      interval_2(interval_2 > (L-1)) = interval_2(interval_2 > (L-1)) - L;
      integer_intervals_2 = [integer_intervals_2, interval_2];

      % satellite 3 위치
      interval_3 = interval_1 + round(L/N*(N-1)); 
      interval_3(interval_3 > (L-1)) = interval_3(interval_3 > (L-1)) - L;
      integer_intervals_3 = [integer_intervals_3, interval_3];

      % satellite 4 위치
      interval_4 = interval_1 + round(L/N*(N-1)); 
      interval_4(interval_4 > (L-1)) = interval_4(interval_4 > (L-1)) - L;
      integer_intervals_4 = [integer_intervals_4, interval_4];

      % satellite 5 위치
      interval_5 = interval_1 + round(L/N*(N-1)); 
      interval_5(interval_5 > (L-1)) = interval_5(interval_5 > (L-1)) - L;
      integer_intervals_5 = [integer_intervals_5, interval_5];

  end


sum_intervals = ...
    [integer_intervals_1, ...
    integer_intervals_2, ...
    integer_intervals_3, ...
    integer_intervals_4, ...
    integer_intervals_5, ...
    % integer_intervals_6, ...
    % integer_intervals_7, ...
    % integer_intervals_8, ...
    % integer_intervals_9, ...
    % integer_intervals_10, ...
    ];

unique_intervals = unique(sort(sum_intervals));

Cov = length(unique_intervals) / L * 100


  if Cov >= 0


        break;
    end
 trials = trials + 1; 
end

% 시뮬레이션 종료 시간 기록
elapsed_time = toc(start_time);

% 결과 출력
fprintf('시행 횟수: %d\n', trials);
fprintf('걸린 시간(초): %.2f\n', elapsed_time);


% 변수 목록
vars = who;

% 'integer_intervals_'로 시작하는 변수들만 필터링
pattern = '^integer_intervals_\d+$'; % 정규 표현식 패턴
N = 0;

for i = 1:length(vars)
    if ~isempty(regexp(vars{i}, pattern, 'once'))
        N = N + 1;
    end
end


%% 그래프

n_1 = [round(L/N * (N-1:-1:0))]+1;

% 동적으로 ff_ 변수를 생성 및 할당
for i = 1:N
    eval(['ff_', num2str(i), ' = ff_0 + n_1(', num2str(i), ');']);
end


% x축 범위 설정
x_range = 0:L-1;

for idx = 1:length(n_1)
    ff = eval(['ff_', num2str(idx)]);
    x_coords1{idx} = []; % 해당 작업공간의 x좌표를 저장할 배열 초기화

    for i = 1:size(ff,1)
        start_time = ff(i, 1);
        end_time = ff(i, 2);

        % y=1일 때의 시작 시간과 종료 시간의 x좌표 계산
        if end_time > L-1 && start_time < L
           start_index = find(x_range >= start_time, 1);
           end_index = find(x_range < end_time-L, 1, 'last');
        elseif start_time > L-1
           start_index = find(x_range >= start_time-L, 1);
           end_index = find(x_range <= end_time-L, 1, 'last');
        else
        start_index = find(x_range >= start_time, 1);
        end_index = find(x_range <= end_time, 1, 'last');
        end

        % 시작 시간부터 종료 시간까지의 x좌표를 배열에 저장
        if end_time > L-1 && start_time < L
          x_coords1{idx} = [x_coords1{idx}, x_range(start_index:L),(0:end_index)];

        else
          x_coords1{idx} = [x_coords1{idx}, x_range(start_index:end_index)];
        end
    end               
end

% 모든 작업공간의 x 좌표를 모은 배열 초기화
all_x_coords1 = [];

%% 각 작업공간의 x 좌표를 all_x_coords1에 추가하고 작은 수부터 오름차순으로 정렬
for idx = 1:length(n_1)
    all_x_coords1 = [all_x_coords1, x_coords1{idx}];
end

all_x_coords1 = sort(all_x_coords1);

% 중복된 값을 제거하고 유일한 값들만 추출
unique_x_coords1 = unique(all_x_coords1);

%%

% 초기화
RV1 = [];
current_value = all_x_coords1(1);
start_index = 1;

% 각 원소에 대해 반복하면서 연속된 값의 구간을 찾음
for i = 2:length(all_x_coords1)
    if all_x_coords1(i) ~= current_value
        % 현재 값과 다른 값이 나타난 경우, 해당 구간의 시작과 끝을 저장하고 현재 값 갱신
        end_index = i - 1;
        RV1 = [RV1; current_value, end_index - start_index + 1];
        current_value = all_x_coords1(i);
        start_index = i;
    end
end

% 마지막 구간에 대한 처리
end_index = length(all_x_coords1);
RV1 = [RV1; current_value, end_index - start_index+1];


%% 빠진 숫자 찾기

% 주어진 행렬
matrix1 = RV1 ;

% 현재 1열에서 빠진 숫자를 찾기
missing_numbers = setdiff(0:L-1, matrix1(:,1));

% 새롭게 추가된 행에 대해 0으로 입력
new_rows = zeros(length(missing_numbers), 2);
new_rows(:, 1) = missing_numbers;

% 새로운 행렬 생성
new_matrix1 = [matrix1; new_rows];

% 1열을 기준으로 행렬을 정렬
sorted_matrix1 = sortrows(new_matrix1, 1);

% Acess Profile

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
    if ff_0(i) < 0
        start_time_new = L + ff_0(i,1);
        end_time_new = L + ff_0(i,2);
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
    if ff_0(i) < 0
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

% x_range = 0:L-1;

subplot(3, 1, 2);

    % y 값 설정
    y_values1p = zeros(size(x_range));
    y_values1p(n_1(1,:)) = 1;

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
plot(sorted_matrix1(:,1), sorted_matrix1(:,2), 'Color', [red1, green1, blue1], 'LineWidth', 1, 'MarkerSize', 1);
title('Coverage Timeline');
hold on;
area(sorted_matrix1(:,1), sorted_matrix1(:,2), 'FaceColor', [red11, green11, blue11]);
plot(sorted_matrix1(:,1), sorted_matrix1(:,2), 'Color', [red1, green1, blue1], 'LineWidth', 1);

xlim([0, L-1]);
ylim([0, 3]);

% 그래프 세부 설정
xlabel('n');
ylabel('bj[n]');
grid on;
xlim([0, L-1]);
