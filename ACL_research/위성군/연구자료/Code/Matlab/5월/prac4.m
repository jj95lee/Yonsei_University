% 초기화
clc;
clear all;

L = 717;
shift = 0;

% 시작점과 끝점 정의
ff_00 = [
    1422.34375	3882.988281;
    12744.29688	15186.19141;
    23907.98828	26012.96875;
    35150.78125	35906.23047;
    55517.20703	56159.45313;
    65349.64844	67430.07813;
    76151.81641	78589.90234
];

ff_0 = floor(ff_00 / (86024 / L)) - shift;

start_points = ff_0(:, 1);
end_points = ff_0(:, 2);

% 결과 초기화
result = [];

% 첫 번째 반복문
for n1 = 0:716
    % 두 번째 반복문
    for n2 = 0:716
        integer_intervals_1 = [];
        integer_intervals_2 = [];
        
        % 구간 생성
        for i = 1:length(ff_0)
            interval_1 = start_points(i) : end_points(i);
            interval_1(interval_1 > 716) = interval_1(interval_1 > 716) - 717;
            integer_intervals_1 = [integer_intervals_1, interval_1];

            interval_2 = interval_1 + n1; 
            interval_2(interval_2 > 716) = interval_2(interval_2 > 716) - 717;
            integer_intervals_2 = [integer_intervals_2, interval_2 + n2];
        end

        % 두 행렬 합치고 중복 제거
        sum_intervals = [integer_intervals_1, integer_intervals_2];
        unique_intervals = unique(sort(sum_intervals));

        % 차이 계산
        difference = length(sum_intervals) - length(unique_intervals);

        % 결과 저장
        result = [result; n1, n2, difference];
    end
end

% 결과에서 중복이 가장 적은 값을 찾기
[min_diff, min_idx] = min(result(:, 3));
optimal_combination = result(min_idx, 1:2);

fprintf('Optimal combination: n1 = %d, n2 = %d with minimum duplicates of %d\n', ...
    optimal_combination(1), optimal_combination(2), min_diff);
