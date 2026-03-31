clc
clear all

trials = 0;
start_time = tic;
L = 720;

while true

% Seoul

ff_00 = ...
[...
36399.19922	37581.73828	;
43423.41797	44669.47266	;
50630.03906	51286.93359	;
70706.30859	71445.05859	;
77369.23828	78631.64063	;
84492.42188	85639.27734
...		
];

% ff_0 -> n 으로 변환
ff_0 = floor(ff_00/(86400/L));

% 구간의 시작점과 끝점을 변수에 저장
start_points = ff_0(:, 1);
end_points = ff_0(:, 2);


%% Random (ex. 위성 22대)

integer_intervals_1 = [];
integer_intervals_2 = [];
integer_intervals_3 = [];
integer_intervals_4 = [];
integer_intervals_5 = [];
integer_intervals_6 = [];
integer_intervals_7 = [];
integer_intervals_8 = [];
integer_intervals_9 = [];
integer_intervals_10 = [];
integer_intervals_11 = [];
integer_intervals_12 = [];
integer_intervals_13 = [];
integer_intervals_14 = [];
integer_intervals_15 = [];
integer_intervals_16 = [];
integer_intervals_17 = [];
integer_intervals_18 = [];
integer_intervals_19 = [];
integer_intervals_20 = [];
integer_intervals_21 = [];
integer_intervals_22 = [];
% integer_intervals_23 = [];
% integer_intervals_24 = [];
% integer_intervals_25 = [];
% integer_intervals_26 = [];
% integer_intervals_27 = [];
% integer_intervals_28 = [];
% integer_intervals_29 = [];
% integer_intervals_30 = [];
% integer_intervals_31 = [];
% integer_intervals_32 = [];
% integer_intervals_33 = [];
% integer_intervals_34 = [];
% integer_intervals_35 = [];
% integer_intervals_36 = [];
% integer_intervals_37 = [];

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
      % satellite 1 위치(Seed Satellite)
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
      interval_4 = interval_1 + round(L/N*(N-2)); 
      interval_4(interval_4 > (L-1)) = interval_4(interval_4 > (L-1)) - L;
      integer_intervals_4 = [integer_intervals_4, interval_4];

      % satellite 5 위치
      interval_5 = interval_1 + round(L/N*(N-3)); 
      interval_5(interval_5 > (L-1)) = interval_5(interval_5 > (L-1)) - L;
      integer_intervals_5 = [integer_intervals_5, interval_5];

      % satellite 6 위치
      interval_6 = interval_1 + round(L/N*(N-4)); 
      interval_6(interval_6 > (L-1)) = interval_6(interval_6 > (L-1)) - L;
      integer_intervals_6 = [integer_intervals_6, interval_6];

      % satellite 7 위치
      interval_7 = interval_1 + round(L/N*(N-5)); 
      interval_7(interval_7 > (L-1)) = interval_7(interval_7 > (L-1)) - L;
      integer_intervals_7 = [integer_intervals_7, interval_7];

      % satellite 8 위치
      interval_8 = interval_1 + round(L/N*(N-6)); 
      interval_8(interval_8 > (L-1)) = interval_8(interval_8 > (L-1)) - L;
      integer_intervals_8 = [integer_intervals_8, interval_8];

      % satellite 9 위치
      interval_9 = interval_1 + round(L/N*(N-7)); 
      interval_9(interval_9 > (L-1)) = interval_9(interval_9 > (L-1)) - L;
      integer_intervals_9 = [integer_intervals_9, interval_9];

      % satellite 10 위치
      interval_10 = interval_1 + round(L/N*(N-8)); 
      interval_10(interval_10 > (L-1)) = interval_10(interval_10 > (L-1)) - L;
      integer_intervals_10 = [integer_intervals_10, interval_10];

      % satellite 11 위치 
      interval_11 = interval_1 + round(L/N*(N-9)); 
      interval_11(interval_11 > (L-1)) = interval_11(interval_11 > (L-1)) - L;
      integer_intervals_11 = [integer_intervals_11, interval_11];

      % satellite 12 위치
      interval_12 = interval_1 + round(L/N*(N-10)); 
      interval_12(interval_12 > (L-1)) = interval_12(interval_12 > (L-1)) - L;
      integer_intervals_12 = [integer_intervals_12, interval_12];

      % satellite 13 위치
      interval_13 = interval_1 + round(L/N*(N-11)); 
      interval_13(interval_13 > (L-1)) = interval_13(interval_13 > (L-1)) - L;
      integer_intervals_13 = [integer_intervals_13, interval_13];

      % satellite 14 위치
      interval_14 = interval_1 + round(L/N*(N-12)); 
      interval_14(interval_14 > (L-1)) = interval_14(interval_14 > (L-1)) - L;
      integer_intervals_14 = [integer_intervals_14, interval_14];

      % satellite 15 위치
      interval_15 = interval_1 + round(L/N*(N-13)); 
      interval_15(interval_15 > (L-1)) = interval_15(interval_15 > (L-1)) - L;
      integer_intervals_15 = [integer_intervals_15, interval_15];

      % satellite 16 위치
      interval_16 = interval_1 + round(L/N*(N-14)); 
      interval_16(interval_16 > (L-1)) = interval_16(interval_16 > (L-1)) - L;
      integer_intervals_16 = [integer_intervals_16, interval_16];

      % satellite 17 위치
      interval_17 = interval_1 + round(L/N*(N-15)); 
      interval_17(interval_17 > (L-1)) = interval_17(interval_17 > (L-1)) - L;
      integer_intervals_17 = [integer_intervals_17, interval_17];

      % satellite 18 위치 
      interval_18 = interval_1 + round(L/N*(N-16)); 
      interval_18(interval_18 > (L-1)) = interval_18(interval_18 > (L-1)) - L;
      integer_intervals_18 = [integer_intervals_18, interval_18];

      % satellite 19 위치
      interval_19 = interval_1 + round(L/N*(N-17)); 
      interval_19(interval_19 > (L-1)) = interval_19(interval_19 > (L-1)) - L;
      integer_intervals_19 = [integer_intervals_19, interval_19];

      % satellite 20 위치 
      interval_20 = interval_1 + round(L/N*(N-18)); 
      interval_20(interval_20 > (L-1)) = interval_20(interval_20 > (L-1)) - L;
      integer_intervals_20 = [integer_intervals_20, interval_20];

      % satellite 21 위치
      interval_21 = interval_1 + round(L/N*(N-19)); 
      interval_21(interval_21 > (L-1)) = interval_21(interval_21 > (L-1)) - L;
      integer_intervals_21 = [integer_intervals_21, interval_21];

      % satellite 22 위치 
      interval_22 = interval_1 + round(L/N*(N-20)); 
      interval_22(interval_22 > (L-1)) = interval_22(interval_22 > (L-1)) - L;
      integer_intervals_22 = [integer_intervals_22, interval_22];

      % % satellite 23 위치
      % interval_23 = interval_1 + round(L/N*(N-21)); 
      % interval_23(interval_23 > (L-1)) = interval_23(interval_23 > (L-1)) - L;
      % integer_intervals_23 = [integer_intervals_23, interval_23];
      % 
      % % satellite 24 위치 
      % interval_24 = interval_1 + round(L/N*(N-22)); 
      % interval_24(interval_24 > (L-1)) = interval_24(interval_24 > (L-1)) - L;
      % integer_intervals_24 = [integer_intervals_24, interval_24];

      % % satellite 25 위치
      % interval_25 = interval_1 + round(L/N*(N-23)); 
      % interval_25(interval_25 > (L-1)) = interval_25(interval_25 > (L-1)) - L;
      % integer_intervals_25 = [integer_intervals_25, interval_25];
      % 
      % % satellite 26 위치
      % interval_26 = interval_1 + round(L/N*(N-24)); 
      % interval_26(interval_26 > (L-1)) = interval_26(interval_26 > (L-1)) - L;
      % integer_intervals_26 = [integer_intervals_26, interval_26];
      % 
      % % satellite 27 위치
      % interval_27 = interval_1 + round(L/N*(N-25)); 
      % interval_27(interval_27 > (L-1)) = interval_27(interval_27 > (L-1)) - L;
      % integer_intervals_27 = [integer_intervals_27, interval_27];
      % 
      % % satellite 28 위치
      % interval_28 = interval_1 + round(L/N*(N-26)); 
      % interval_28(interval_28 > (L-1)) = interval_28(interval_28 > (L-1)) - L;
      % integer_intervals_28 = [integer_intervals_28, interval_28];
      % 
      % % satellite 29 위치
      % interval_29 = interval_1 + round(L/N*(N-27)); 
      % interval_29(interval_29 > (L-1)) = interval_29(interval_29 > (L-1)) - L;
      % integer_intervals_29 = [integer_intervals_29, interval_29];

      % % satellite 30 위치 
      % interval_30 = interval_1 + round(L/N*(N-28)); 
      % interval_30(interval_30 > (L-1)) = interval_30(interval_30 > (L-1)) - L;
      % integer_intervals_30 = [integer_intervals_30, interval_30];
  end


sum_intervals = ...
    [integer_intervals_1, ...
    integer_intervals_2, ...
    integer_intervals_3, ...
    integer_intervals_4, ...
    integer_intervals_5, ...
    integer_intervals_6, ...
    integer_intervals_7, ...
    integer_intervals_8, ...
    integer_intervals_9, ...
    integer_intervals_10, ...
    integer_intervals_11, ...
    integer_intervals_12, ...
    integer_intervals_13, ...
    integer_intervals_14, ...
    integer_intervals_15, ...
    integer_intervals_16, ...
    integer_intervals_17, ...
    integer_intervals_18, ...
    integer_intervals_19, ...
    integer_intervals_20, ...
    integer_intervals_21, ...
    integer_intervals_22, ...
    % integer_intervals_23, ...
    % integer_intervals_24, ...
    % integer_intervals_25, ...
    % integer_intervals_26, ...
    % integer_intervals_27, ...
    % integer_intervals_28, ...
    % integer_intervals_29, ...
    % integer_intervals_30, ...
    ];

unique_intervals = unique(sort(sum_intervals));

% 최소 및 최대 값 계산
min_val = min(sum_intervals);
max_val = max(sum_intervals);

% 모든 숫자에 대해 개수를 기록할 행렬 초기화
sum_intervals_ac = zeros(L, 2);

% 1열에 숫자들 채우기
sum_intervals_ac(:, 1) = (0:L-1)';

% 각 숫자의 개수 세기
for i = min_val:max_val
    sum_intervals_ac(i - min_val + 1, 2) = sum(sum_intervals == i);
end


Cov = length(unique_intervals) / L * 100

  if Cov >= 50

        break;
    end
 trials = trials + 1; 
end


% 시뮬레이션 종료 시간 기록
elapsed_time = toc(start_time);

% 결과 출력
fprintf('시행 횟수: %d\n', trials);
fprintf('걸린 시간(초): %.2f\n', elapsed_time);


%% 그래프

n_1 = [round(L/N * (0:1:N-1))]+1;


% ff_ 변수를 생성 및 할당
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

        % 시작 시간과 종료 시간의 x좌표 계산(L의 위치에 따른 보정)
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

%% Acess Profile

% x축 범위 설정
x_range = 0:L-1;

figure;
subplot(3, 1, 1);

hex_color1 = '#76A646';

red1 = hex2dec(hex_color1(2:3)) / 255; % R 값
green1 = hex2dec(hex_color1(4:5)) / 255; % G 값
blue1 = hex2dec(hex_color1(6:7)) / 255; % B 값

hex_color11 = '#76A646';

red11 = hex2dec(hex_color11(2:3)) / 255; % R 값
green11 = hex2dec(hex_color11(4:5)) / 255; % G 값
blue11 = hex2dec(hex_color11(6:7)) / 255; % B 값


% 각 구간별로 y=1인 위치 표시
for i = 1:size(ff_0, 1)

    start_time = ff_0(i, 1);
    end_time = ff_0(i, 2);

    % 시작 시간과 종료 시간에 해당하는 x 인덱스 찾기
    start_index = find(x_range >= start_time, 1);
    end_index = find(x_range <= end_time, 1, 'last');

    % y 값 설정
    y_values1a = zeros(size(x_range));
    y_values1a(start_index:end_index) = 1;

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
yticks([0 1]);

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
yticks([0 1]);

% 그래프 세부 설정
xlabel('n');
ylabel('x[n]');
grid on;
xlim([0, L-1]);
ylim([0, 1]);
yticks([0 1]);

%% Coverage Timline

x_range = 0:L-1;

subplot(3, 1, 3);
plot(sum_intervals_ac(:,1), sum_intervals_ac(:,2), 'Color', [red1, green1, blue1], 'LineWidth', 1, 'MarkerSize', 1);
title('Coverage Timeline');
hold on;
area(sum_intervals_ac(:,1), sum_intervals_ac(:,2), 'FaceColor', [red11, green11, blue11]);
plot(sum_intervals_ac(:,1), sum_intervals_ac(:,2), 'Color', [red1, green1, blue1], 'LineWidth', 1);

xlim([0, L-1]);
ylim([0, 3]);

% 그래프 세부 설정
xlabel('n');
ylabel('bj[n]');
grid on;
xlim([0, L-1]);
