clc
clear all


% 주어진 행렬 ff_0
ff_0 = [235	261
362	397
499	534
635	660];



% 구간의 시작점과 끝점을 변수에 저장
start_points = ff_0(:, 1);
end_points = ff_0(:, 2);

% 모든 구간을 정수로 표시해서 작업 공간에 저장
result = [];

for n = 0:716
integer_intervals_1 = [];
integer_intervals_2 = [];
integer_intervals_3 = [];
integer_intervals_4 = [];
integer_intervals_5 = [];
integer_intervals_6 = [];
integer_intervals_7 = [];
integer_intervals_8 = [];




%%    519   617   419   192   228   180   109

for i = 1:length(ff_0)
    interval_1 = start_points(i):end_points(i) ;
    interval_1(interval_1 > 716) = interval_1(interval_1 > 716) - 717; 
    integer_intervals_1 = [integer_intervals_1, interval_1];

    % interval_2 = interval_1 + 519; 
    % interval_2(interval_2 > 716) = interval_2(interval_2 > 716) - 717;
    % integer_intervals_2 = [integer_intervals_2, interval_2];
    % 
    % interval_3 = interval_1 + n; 
    % interval_3(interval_3 > 716) = interval_3(interval_3 > 716) - 717;
    % integer_intervals_3 = [integer_intervals_3, interval_3];

    % interval_4 = interval_1 + n; 
    % interval_4(interval_4 > 716) = interval_4(interval_4 > 716) - 717;
    % integer_intervals_4 = [integer_intervals_4, interval_4];
    % 
    % interval_5 = interval_1 + n; 
    % interval_5(interval_5 > 716) = interval_5(interval_5 > 716) - 717;
    % integer_intervals_5 = [integer_intervals_5, interval_5];
    % 
    % interval_6 = interval_1 + n; 
    % interval_6(interval_6 > 716) = interval_6(interval_6 > 716) - 717;
    % integer_intervals_6 = [integer_intervals_6, interval_6];
    % 
    % interval_7 = interval_1 + n; 
    % interval_7(interval_7 > 716) = interval_7(interval_7 > 716) - 717;
    % integer_intervals_7 = [integer_intervals_7, interval_7];
    % 
    % interval_8 = interval_1 + n; 
    % interval_8(interval_8 > 716) = interval_8(interval_8 > 716) - 717;
    % integer_intervals_8 = [integer_intervals_8, interval_8];


end




%%


% 두 행렬을 합쳐서 하나의 행렬로 만듭니다.
sum_intervals = ...
    [integer_intervals_1, ...
    integer_intervals_2, ...
    integer_intervals_3, ...
    integer_intervals_4, ...
    integer_intervals_5, ...
    integer_intervals_6, ...
    integer_intervals_7, ...
    ];

% 행렬을 오름차순으로 정렬하고 중복된 요소를 제거합니다.
unique_intervals = unique(sort(sum_intervals));


% 두 행렬의 길이 차이 계산
difference = length(sum_intervals) - length(unique_intervals) ;

 % 결과를 작업 공간에 저장
    result = [result; n, difference];
end


A = length(unique_intervals) / 717 * 100

% result 변수에서 2열이 최소인 행을 찾습니다.
rows_with_min = result(result(:, 2) == min(result(:, 2)), :);

k_1 = rows_with_min(:,1);
