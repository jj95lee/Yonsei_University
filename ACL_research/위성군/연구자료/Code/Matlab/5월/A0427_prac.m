clear all

% 주어진 ff_0
ff_00 = ...
[...
1422.34375	3882.988281
12744.29688	15186.19141
23907.98828	26012.96875
35150.78125	35906.23047
55517.20703	56159.45313
65349.64844	67430.07813
76151.81641	78589.90234
...		
];

L = 717;
shift = 0;
ff_0 = floor(ff_00 / (86024 / L)) - shift;
start_points = ff_0(:, 1);
end_points = ff_0(:, 2);

% 첫 번째 랜덤 시도
result_1 = [];

for n = 0:716
    integer_intervals_1 = [];
    integer_intervals_2 = [];
    
    for i = 1:length(ff_0)
        interval_1 = start_points(i) : end_points(i);
        interval_1(interval_1 > 716) = interval_1(interval_1 > 716) - 717;
        integer_intervals_1 = [integer_intervals_1, interval_1];
        
        interval_2 = interval_1 + n;
        interval_2(interval_2 > 716) = interval_2(interval_2 > 716) - 717;
        integer_intervals_2 = [integer_intervals_2, interval_2];
    end
    
    sum_intervals = [integer_intervals_1, integer_intervals_2];
    unique_intervals = unique(sort(sum_intervals));
    difference = length(sum_intervals) - length(unique_intervals);
    result_1 = [result_1; n, difference];
end

% k_1 구하기
min_diff = min(result_1(:, 2));
k_1 = result_1(result_1(:, 2) == min_diff, 1);

% k_1에 대해 두 번째 랜덤 시도
k_2_results = {}; % 모든 k_2 결과를 저장하기 위한 셀 배열

for k_val = 1:length(k_1)
    n = k_1(k_val); % 첫 번째 랜덤 시도의 k_1 값
    
    result_2 = [];
    
    for j = 0:716
        integer_intervals_1 = [];
        integer_intervals_2 = [];
        integer_intervals_3 = [];
        
        for i = 1:length(ff_0)
            interval_1 = start_points(i) : end_points(i);
            interval_1(interval_1 > 716) = interval_1(interval_1 > 716) - 717;
            integer_intervals_1 = [integer_intervals_1, interval_1];
            
            interval_2 = interval_1 + n; % n은 첫 번째 시행에서의 k_1 값
            interval_2(interval_2 > 716) = interval_2(interval_2 > 716) - 717;
            integer_intervals_2 = [integer_intervals_2, interval_2];
            
            interval_3 = interval_1 + j;
            interval_3(interval_3 > 716) = interval_3(interval_3 > 716) - 717;
            integer_intervals_3 = [integer_intervals_3, interval_3];
        end
        
        sum_intervals = [integer_intervals_1, integer_intervals_2, integer_intervals_3];
        unique_intervals = unique(sort(sum_intervals));
        difference = length(sum_intervals) - length(unique_intervals);
        result_2 = [result_2; j, difference];
    end
    
    min_diff_2 = min(result_2(:, 2)); % 두 번째 시행의 최소 중복값
    k_2 = result_2(result_2(:, 2) == min_diff_2, 1); % 최소 중복을 만드는 k_2 값들
    k_2_results{k_val} = k_2; % 결과를 셀 배열에 저장
end

% 모든 k_2 행렬 출력
for k_val = 1:length(k_2_results)
    fprintf('k_1 값 %d에 대한 k_2 행렬:\n', k_1(k_val));
    disp(k_2_results{k_val}); % 각 k_1에 대한 k_2 결과를 출력
end

% 셀 배열에서 각 행의 최대 길이 찾기
max_length = 0;

for k_val = 1:length(k_2_results)
    current_length = length(k_2_results{k_val});
    if current_length > max_length
        max_length = current_length;
    end
end

% 144*320 행렬 생성 (320은 k_2_results의 열 수)
k_2_matrix = zeros(max_length, length(k_2_results));

% 셀 배열에서 행렬로 데이터 복사
for k_val = 1:length(k_2_results)
    k_2_values = k_2_results{k_val};
    k_2_matrix(1:length(k_2_values), k_val) = k_2_values;
end


