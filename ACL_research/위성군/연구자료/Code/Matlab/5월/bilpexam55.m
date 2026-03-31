clear all

L = 717;

%% 저궤도 

% n의 배열
n_1 = [65, 144, 285, 361] + 1 ;

NumberofSatellite = 4;
ff_00 = ...
[...
1422.34375	3882.988281
12744.29688	15186.19141
23907.98828	26012.96875
35150.78125	35906.23047
55517.20703	56159.45313
65349.64844	67430.07813
76151.81641	78589.90234...		
];

% ff_0= [...
% 11+L 32+L
% 106+L 126+L
% 199+L 216+L
% 292+L 299+L
% 462	467
% 544	561
% 634	654]-355;

% ff_0 -> n 으로 변환
ff_0 = floor(ff_00/(86024/717))-358;

% 음수인 요소 찾기
negative_indices = ff_0 < 0;

ff_0(negative_indices) = ff_0(negative_indices) + L;

%ff_0 = seedsatellite
ff_1 = ff_0 + n_1(1);
ff_2 = ff_0 + n_1(2);
ff_3 = ff_0 + n_1(3);
ff_4 = ff_0 + n_1(4);
% ff_5 = ff_0 + n_1(5);
% ff_6 = ff_0 + n_1(6);
% ff_7 = ff_0 + n_1(7);
% ff_8 = ff_0 + n_1(8);
% ff_9 = ff_0 + n_1(9);
% ff_10 = ff_0 + n_1(10);


% x축의 범위 설정 (0부터 716까지)
x_range =  0:1:716;

%% 모든 작업공간 그래프 그리기
% for idx = 1:NumberofSatellite
%     ff = eval(['ff_', num2str(idx)]);
% 
%     for i = 1:size(ff,1)
%         start_time = ff(i, 1);
%         end_time = ff(i, 2);
% 
%         start_index = find(x_range >= start_time, 1);
%         end_index = find(x_range <= end_time, 1, 'last');
% 
%         y_range = zeros(size(x_range));
%         y_range(start_index:end_index) = 1;
% 
% 
% % 첫 번째 구간 추가
% range1 = x_range(start_index:end_index);
% 
% % 두 번째 구간 추가 (만약 end_time이 717을 초과하는 경우)
% if end_time > 716 & start_time < 717
%     range2 = [x_range(1:end_time-716)];
% elseif start_time > 716
%     range2 = [x_range(start_time-716:end_time-716)];
%  else
%     range2 = range1;
% end
% 
% 
% % 그래프 그리기
%  fill([x_range(start_index:end_index), fliplr(x_range(start_index:end_index))], ...
%              [y_range(start_index:end_index), zeros(1, length(start_index:end_index))], ...
%              [0.3, 0.6, 0.9], 'EdgeColor', 'none');
%         hold on;
% 
%  fill([range2, fliplr(range2)], ...
%              [ones(1, length(range2)), zeros(1, length(range2))], ...
%              [0.3, 0.6, 0.9], 'EdgeColor', 'none');
%     end
% end
% 
% % 그래프 세부 설정
% xlabel('시간 (초)');
% ylabel('Coverage');
% title('구간 그래프');
% grid on;
% xlim([0, 716]);

%%

for idx = 1:NumberofSatellite
    ff = eval(['ff_', num2str(idx)]);
    x_coords1{idx} = []; % 해당 작업공간의 x좌표를 저장할 배열 초기화

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
          x_coords1{idx} = [x_coords1{idx}, x_range(start_index:717),(0:end_index)];
      
        else
          x_coords1{idx} = [x_coords1{idx}, x_range(start_index:end_index)];
        end
        


    end               
end

% 모든 작업공간의 x 좌표를 모은 배열 초기화
all_x_coords1 = [];

%% 각 작업공간의 x 좌표를 all_x_coords1에 추가하고 작은 수부터 오름차순으로 정렬
for idx = 1:NumberofSatellite
    all_x_coords1 = [all_x_coords1, x_coords1{idx}];
end

all_x_coords1 = sort(all_x_coords1);

% 중복된 값을 제거하고 유일한 값들만 추출
unique_x_coords1 = unique(all_x_coords1);

coverage1 = length(unique_x_coords1) / 717 * 100;

%%

% 초기화
RV = [];
current_value = all_x_coords1(1);
start_index = 1;

% 각 원소에 대해 반복하면서 연속된 값의 구간을 찾음
for i = 2:length(all_x_coords1)
    if all_x_coords1(i) ~= current_value
        % 현재 값과 다른 값이 나타난 경우, 해당 구간의 시작과 끝을 저장하고 현재 값 갱신
        end_index = i - 1;
        RV = [RV; current_value, end_index - start_index + 1];
        current_value = all_x_coords1(i);
        start_index = i;
    end
end

% 마지막 구간에 대한 처리
end_index = length(all_x_coords1);
RV = [RV; current_value, end_index - start_index + 1];


%% 빠진 숫자 찾기

% 주어진 행렬
matrix = RV ;

% 현재 1열에서 빠진 숫자를 찾기
missing_numbers = setdiff(0:716, matrix(:,1));

% 새롭게 추가된 행에 대해 0으로 입력
new_rows = zeros(length(missing_numbers), 2);
new_rows(:, 1) = missing_numbers;

% 새로운 행렬 생성
new_matrix = [matrix; new_rows];

% 1열을 기준으로 행렬을 정렬
sorted_matrix1 = sortrows(new_matrix, 1);


% %% Acess Profile
% % x축 범위 설정
% x_range = 0:716;
% 
% figure;
% subplot(3, 1, 1);
% 
hex_color1 = '#005eb5';

r1 = hex2dec(hex_color1(2:3)) / 255; % R 값
g1 = hex2dec(hex_color1(4:5)) / 255; % G 값
b1 = hex2dec(hex_color1(6:7)) / 255; % B 값

hex_color11 = '#6ca6d6';

r11 = hex2dec(hex_color11(2:3)) / 255; % R 값
g11 = hex2dec(hex_color11(4:5)) / 255; % G 값
b11 = hex2dec(hex_color11(6:7)) / 255; % B 값
% 
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
    % y 값 설정
    y_values1a = zeros(size(x_range));
    y_values1a(start_index:end_index) = 1;
% 
%     % 그래프 그리기
%     plot(x_range, y_values1a, '-o r', 'LineWidth', 1, 'MarkerSize', 1, 'Color', [r1, g1, b1]);
%     hold on;
%     area(x_range, y_values1a, 'FaceColor', [r11, g11, b11]);
%     plot(x_range, y_values1a, 'Color', [r1, g1, b1], 'LineWidth', 1.5); 
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
    % y 값 설정
    y_values1p = zeros(size(x_range));
    y_values1p(n_1(1,:)) = 1;
% 
% 
% plot(x_range, y_values1p, 'Color', [r1, g1, b1], 'LineWidth', 1, 'MarkerSize', 1);
% title('Constellation Pattern Vector');
% hold on;
% 
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
% plot(sorted_matrix1(:,1), sorted_matrix1(:,2), 'Color', [r1, g1, b1], 'LineWidth', 1, 'MarkerSize', 1);
% title('Coverage Timeline');
% hold on;
% area(sorted_matrix1(:,1), sorted_matrix1(:,2), 'FaceColor', [r11, g11, b11]);
% plot(sorted_matrix1(:,1), sorted_matrix1(:,2), 'Color', [r1, g1, b1], 'LineWidth', 1);
% 
% xlim([0, 716]);
% ylim([0, 4]);
% 
% % 그래프 세부 설정
% xlabel('n');
% ylabel('bj[n]');
% grid on;
% xlim([0, 716]);


%% V
% 작업공간을 저장할 변수 초기화
acp = zeros(1, 717);

% 각 구간에 대해 1로 설정
for i = 1:size(ff_0, 1)
    start_index = ff_0(i, 1);
    end_index = ff_0(i, 2);
    acp(start_index:end_index) = 1;
end


shiftAmount = 1;

V = [];

for i = 1:length(acp)
    newRow = circshift(acp, [0, i-shiftAmount]);
    V = [V; newRow];
end


%% 고궤도

% n의 배열
n_2 = [208 428 523 608 634 702]+1 ;

L = 717;

NumberofSatellite = 6;

gg_00 = ...
[...
1720.332031	5496.875
17319.04297	20498.98438
72635.99609	75944.84375...		
];


% gg_0 -> n 으로 변환
gg_0 = floor(gg_00/(86024/717))-358;

% 음수인 요소 찾기
negative_indices = gg_0 < 0;

gg_0(negative_indices) = gg_0(negative_indices) + L;

%gg_0 = seedsatellite
gg_1 = gg_0 + n_2(1);
gg_2 = gg_0 + n_2(2);
gg_3 = gg_0 + n_2(3);
gg_4 = gg_0 + n_2(4);
gg_5 = gg_0 + n_2(5);
gg_6 = gg_0 + n_2(6);
% gg_7 = gg_0 + n_2(7);
% gg_8 = gg_0 + n_2(8);
% gg_9 = gg_0 + n_2(9);
% gg_10 = gg_0 + n_2(10);


% 시작 시각과 종료 시각이 저장된 행렬
start_time = gg_0(:, 1);
end_time = gg_0(:, 2);

%%

for idx = 1:NumberofSatellite
    gg = eval(['gg_', num2str(idx)]);
    x_coords2{idx} = []; % 해당 작업공간의 x좌표를 저장할 배열 초기화

    for i = 1:size(gg,1)
        start_time = gg(i, 1);
        end_time = gg(i, 2);

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
          x_coords2{idx} = [x_coords2{idx}, x_range(start_index:717),(0:end_index)];
      
        else
          x_coords2{idx} = [x_coords2{idx}, x_range(start_index:end_index)];
        end
        


    end               
end

% 모든 작업공간의 x 좌표를 모은 배열 초기화
all_x_coords2 = [];

%%

for idx = 1:NumberofSatellite
    gg = eval(['gg_', num2str(idx)]);
    x_coords2{idx} = []; % 해당 작업공간의 x좌표를 저장할 배열 초기화

    for i = 1:size(gg,1)
        start_time = gg(i, 1);
        end_time = gg(i, 2);

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
          x_coords2{idx} = [x_coords2{idx}, x_range(start_index:717),(0:end_index)];
      
        else
          x_coords2{idx} = [x_coords2{idx}, x_range(start_index:end_index)];
        end
        


    end               
end

% 모든 작업공간의 x 좌표를 모은 배열 초기화
all_x_coords2 = [];

%% 각 작업공간의 x 좌표를 all_x_coords2에 추가하고 작은 수부터 오름차순으로 정렬
for idx = 1:NumberofSatellite
    all_x_coords2 = [all_x_coords2, x_coords2{idx}];
end

all_x_coords2 = sort(all_x_coords2);

% 중복된 값을 제거하고 유일한 값들만 추출
unique_x_coords2 = unique(all_x_coords2);

coverage2 = length(unique_x_coords2) / 717 * 100;

%%

% 초기화
RV = [];
current_value = all_x_coords2(1);
start_index = 1;

% 각 원소에 대해 반복하면서 연속된 값의 구간을 찾음
for i = 2:length(all_x_coords2)
    if all_x_coords2(i) ~= current_value
        % 현재 값과 다른 값이 나타난 경우, 해당 구간의 시작과 끝을 저장하고 현재 값 갱신
        end_index = i - 1;
        RV = [RV; current_value, end_index - start_index + 1];
        current_value = all_x_coords2(i);
        start_index = i;
    end
end

% 마지막 구간에 대한 처리
end_index = length(all_x_coords2);
RV = [RV; current_value, end_index - start_index + 1];


%% 빠진 숫자 찾기

% 주어진 행렬
matrix = RV ;

% 현재 1열에서 빠진 숫자를 찾기
missing_numbers = setdiff(0:716, matrix(:,1));

% 새롭게 추가된 행에 대해 0으로 입력
new_rows = zeros(length(missing_numbers), 2);
new_rows(:, 1) = missing_numbers;

% 새로운 행렬 생성
new_matrix = [matrix; new_rows];

% 1열을 기준으로 행렬을 정렬
sorted_matrix2 = sortrows(new_matrix, 1);


% %% Acess Profile
% % x축 범위 설정
% x_range = 0:716;
% 
% figure;
% subplot(3, 1, 1);
% 
hex_color2 = '#cc4e01';

r2 = hex2dec(hex_color2(2:3)) / 255; % R 값
g2 = hex2dec(hex_color2(4:5)) / 255; % G 값
b2 = hex2dec(hex_color2(6:7)) / 255; % B 값

hex_color22 = '#fe9882';

r22 = hex2dec(hex_color22(2:3)) / 255; % R 값
g22 = hex2dec(hex_color22(4:5)) / 255; % G 값
b22 = hex2dec(hex_color22(6:7)) / 255; % B 값
% 
% 
% % 각 구간별로 y=1인 위치 표시
% for i = 1:size(gg_0, 1)
%     start_time = gg_0(i, 1);
%     end_time = gg_0(i, 2);
% 
%     % 시작 시간과 종료 시간에 해당하는 x 인덱스 찾기
%     start_index = find(x_range >= start_time, 1);
%     end_index = find(x_range <= end_time, 1, 'last');
% 
    % y 값 설정
    y_values2a = zeros(size(x_range));
    y_values2a(start_index:end_index) = 1;
% 
%     % 그래프 그리기
%     plot(x_range, y_values2a, '-o r', 'LineWidth', 1, 'MarkerSize', 1, 'Color', [r2, g2, b2]);
%     hold on;
%     area(x_range, y_values2a, 'FaceColor', [r22, g22, b22]);
%     plot(x_range, y_values2a, 'Color', [r2, g2, b2], 'LineWidth', 1.5);
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
    % y 값 설정
    y_values2p = zeros(size(x_range));
    y_values2p(n_2(1,:)) = 1;
% 
% 
% plot(x_range, y_values2p, 'Color', [r2, g2, b2], 'LineWidth', 1, 'MarkerSize', 1);
% title('Constellation Pattern Vector');
% xlim([0, 716]);
% ylim([0, 1]);
% 
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
% plot(sorted_matrix2(:,1), sorted_matrix2(:,2), 'Color', [r2, g2, b2], 'LineWidth', 1, 'MarkerSize', 1);
% title('Coverage Timeline');
% hold on;
%     area(sorted_matrix2(:,1), sorted_matrix2(:,2), 'FaceColor', [r22, g22, b22]);
%     plot(sorted_matrix2(:,1), sorted_matrix2(:,2), 'Color', [r2, g2, b2], 'LineWidth', 1);
% 
% xlim([0, 716]);
% ylim([0, 4]);
% 
% % 그래프 세부 설정
% xlabel('n');
% ylabel('bj[n]');
% grid on;
% xlim([0, 716]);

%%

% 새로운 그래프에 두 개의 플롯을 함께 보기
figure3 = figure;
subplot(3, 1, 1);

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
    plot(x_range, y_values1a, '-o r', 'LineWidth', 1, 'MarkerSize', 1, 'Color', [r1, g1, b1]);
    hold on;
    area(x_range, y_values1a, 'FaceColor', [r11, g11, b11], 'FaceAlpha', 0.7);
    plot(x_range, y_values1a, 'Color', [r1, g1, b1], 'LineWidth', 1.5);
end

for i = 1:size(gg_0, 1)
    start_time = gg_0(i, 1);
    end_time = gg_0(i, 2);
    
    % 시작 시간과 종료 시간에 해당하는 x 인덱스 찾기
    start_index = find(x_range >= start_time, 1);
    end_index = find(x_range <= end_time, 1, 'last');
    
    % y 값 설정
    y_values2a = zeros(size(x_range));
    y_values2a(start_index:end_index) = 1;
    
    % 그래프 그리기
    plot(x_range, y_values2a, '-o r', 'LineWidth', 1, 'MarkerSize', 1, 'Color', [r2, g2, b2]);
    hold on;
    area(x_range, y_values2a, 'FaceColor', [r22, g22, b22], 'FaceAlpha', 0.7);
    plot(x_range, y_values2a, 'Color', [r2, g2, b2], 'LineWidth', 1.5);
end

hold on; % 여러 플롯을 함께 그리기 위해 hold on 사용
plot(x_range, y_values1a, 'Color', [r1, g1, b1]); % 첫 번째 플롯 복사
plot(x_range, y_values2a, 'Color', [r2, g2, b2]); % 두 번째 플롯 복사
xlim([0, 716]);
ylim([0, 1]);
yticks([0, 1]);

% 그래프 세부 설정
ylabel('vo,j[n]');
title('Seed Satellite Access Profile');
grid on;



subplot(3, 1, 2);
hold on; % 여러 플롯을 함께 그리기 위해 hold on 사용
plot(x_range, y_values1p, 'Color', [r1, g1, b1]); % 첫 번째 플롯 복사
plot(x_range, y_values2p, 'Color', [r2, g2, b2]); % 두 번째 플롯 복사
xlim([0, 716]);
ylim([0, 1]);
yticks([0, 1]);

% 그래프 세부 설정
ylabel('x[n]');
title('Constellation Pattern Vector');
grid on;



subplot(3, 1, 3);
hold on; % 여러 플롯을 함께 그리기 위해 hold on 사용
plot(sorted_matrix1(:,1), sorted_matrix1(:,2), 'Color', [r1, g1, b1]); % 첫 번째 플롯 복사
hold on;
area(sorted_matrix1(:,1), sorted_matrix1(:,2), 'FaceColor', [r11, g11, b11], 'FaceAlpha', 0.7);
plot(sorted_matrix1(:,1), sorted_matrix1(:,2), 'Color', [r1, g1, b1], 'LineWidth', 1.5);

plot(sorted_matrix2(:,1), sorted_matrix2(:,2), 'Color', [r2, g2, b2]); % 두 번째 플롯 복사
hold on;
area(sorted_matrix2(:,1), sorted_matrix2(:,2), 'FaceColor', [r22, g22, b22], 'FaceAlpha', 0.7);
plot(sorted_matrix2(:,1), sorted_matrix2(:,2), 'Color', [r2, g2, b2], 'LineWidth', 1.5);

xlim([0, 716]);
ylim([0, 3]);

% 그래프 세부 설정
ylabel('n');
ylabel('bj[n]');
title('Coverage Timeline');
grid on;

hold off; % hold 해제

sorted_matrix3(:,2) = sorted_matrix1(:,2) + sorted_matrix2(:,2);
count = sum(sorted_matrix3(:, 2) >= 1); 

coverage_total = count / 717 * 100
