clc
clear all

A_18 = 0;
trials = 0;
start_time = tic;

while true


% % 주어진 행렬 ff_0
% ff_0 = [303	313
%         362	372
%         422	428
%         590	595
%         645	655
%         704	714];

% 주어진 행렬 ff_0
ff_0 = [303 313
362	372
422	427
589	595
645	655
704	714];

length_ap = ff_0(1,2)-ff_0(1,1)+ff_0(2,2)-ff_0(2,1)+ff_0(3,2)-ff_0(3,1)+...
+ff_0(4,2)-ff_0(4,1)+ff_0(5,2)-ff_0(5,1)+ff_0(6,2)-ff_0(6,1)+6;



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
result_11 = [];
result_12 = [];
result_13 = [];
result_14 = [];
result_15 = [];
result_16 = [];
result_17 = [];
result_18 = [];
result_19 = [];




%% Random 1

for n = 0:719
integer_intervals_1 = [];
integer_intervals_2 = [];


for i = 1:length(ff_0)
    interval_1 = start_points(i) :end_points(i) ;
    interval_1(interval_1 > 719) = interval_1(interval_1 > 719) - 720; 
    integer_intervals_1 = [integer_intervals_1, interval_1];

    interval_2 = interval_1 + n; 
    interval_2(interval_2 > 719) = interval_2(interval_2 > 719) - 720;
    integer_intervals_2 = [integer_intervals_2, interval_2];

   
end


% 두 행렬을 합쳐서 하나의 행렬로 만듭니다.
sum_intervals = ...
    [integer_intervals_1, ...
    integer_intervals_2, ...
    ];

% 행렬을 오름차순으로 정렬하고 중복된 요소를 제거합니다.
unique_intervals = unique(sort(sum_intervals));


% 두 행렬의 길이 차이 계산
difference = length(sum_intervals) - length(unique_intervals) ;

 % 결과를 작업 공간에 저장
    result_1 = [result_1; n, difference];
end

A_1 = length(unique_intervals) / 720 * 100;

% result 변수에서 2열이 최소인 행을 찾습니다.
rows_with_min_1 = result_1(result_1(:, 2) == min(result_1(:, 2)), :);

k_1 = rows_with_min_1(:,1);


% 랜덤 인덱스 선택
random_index_1 = randi([1, length(k_1)]);

% 선택된 숫자
random_number_1 = k_1(random_index_1);



%% Random 2

for n = 0:719
integer_intervals_1 = [];
integer_intervals_2 = [];
integer_intervals_3 = [];


for i = 1:length(ff_0)
    interval_1 = start_points(i) :end_points(i);
    interval_1(interval_1 > 719) = interval_1(interval_1 > 719) - 720; 
    integer_intervals_1 = [integer_intervals_1, interval_1];

    interval_2 = interval_1 + random_number_1; 
    interval_2(interval_2 > 719) = interval_2(interval_2 > 719) - 720;
    integer_intervals_2 = [integer_intervals_2, interval_2];

    interval_3 = interval_1 + n; 
    interval_3(interval_3 > 719) = interval_3(interval_3 > 719) - 720;
    integer_intervals_3 = [integer_intervals_3, interval_3];

   
end


% 두 행렬을 합쳐서 하나의 행렬로 만듭니다.
sum_intervals = ...
    [integer_intervals_1, ...
    integer_intervals_2, ...
    integer_intervals_3, ...
    ];

% 행렬을 오름차순으로 정렬하고 중복된 요소를 제거합니다.
unique_intervals = unique(sort(sum_intervals));


% 두 행렬의 길이 차이 계산
difference = length(sum_intervals) - length(unique_intervals) ;

 % 결과를 작업 공간에 저장
    result_2 = [result_2; n, difference];
end

A_2 = length(unique_intervals) / 720 * 100;

% result 변수에서 2열이 최소인 행을 찾습니다.
rows_with_min_2 = result_2(result_2(:, 2) == min(result_2(:, 2)), :);

k_2 = rows_with_min_2(:,1);


% 랜덤 인덱스 선택
random_index_2 = randi([1, length(k_2)]);

% 선택된 숫자
random_number_2 = k_2(random_index_2);


%% Random 3

for n = 0:719
integer_intervals_1 = [];
integer_intervals_2 = [];
integer_intervals_3 = [];
integer_intervals_4 = [];


for i = 1:length(ff_0)
    interval_1 = start_points(i):end_points(i);
    interval_1(interval_1 > 719) = interval_1(interval_1 > 719) - 720; 
    integer_intervals_1 = [integer_intervals_1, interval_1];

    interval_2 = interval_1 + random_number_1; 
    interval_2(interval_2 > 719) = interval_2(interval_2 > 719) - 720;
    integer_intervals_2 = [integer_intervals_2, interval_2];

    interval_3 = interval_1 + random_number_2; 
    interval_3(interval_3 > 719) = interval_3(interval_3 > 719) - 720;
    integer_intervals_3 = [integer_intervals_3, interval_3];

    interval_4 = interval_1 + n; 
    interval_4(interval_4 > 719) = interval_4(interval_4 > 719) - 720;
    integer_intervals_4 = [integer_intervals_4, interval_4];

   
end


% 두 행렬을 합쳐서 하나의 행렬로 만듭니다.
sum_intervals = ...
    [integer_intervals_1, ...
    integer_intervals_2, ...
    integer_intervals_3, ...
    integer_intervals_4, ...
    ];

% 행렬을 오름차순으로 정렬하고 중복된 요소를 제거합니다.
unique_intervals = unique(sort(sum_intervals));


% 두 행렬의 길이 차이 계산
difference = length(sum_intervals) - length(unique_intervals) ;

 % 결과를 작업 공간에 저장
    result_3 = [result_3; n, difference];
end

A_3 = length(unique_intervals) / 720 * 100;

% result 변수에서 2열이 최소인 행을 찾습니다.
rows_with_min_3 = result_3(result_3(:, 2) == min(result_3(:, 2)), :);

k_3 = rows_with_min_3(:,1);


% 랜덤 인덱스 선택
random_index_3 = randi([1, length(k_3)]);

% 선택된 숫자
random_number_3 = k_3(random_index_3);

%% Random 4

for n = 0:719
integer_intervals_1 = [];
integer_intervals_2 = [];
integer_intervals_3 = [];
integer_intervals_4 = [];
integer_intervals_5 = [];


for i = 1:length(ff_0)
    interval_1 = start_points(i):end_points(i);
    interval_1(interval_1 > 719) = interval_1(interval_1 > 719) - 720; 
    integer_intervals_1 = [integer_intervals_1, interval_1];

    interval_2 = interval_1 + random_number_1; 
    interval_2(interval_2 > 719) = interval_2(interval_2 > 719) - 720;
    integer_intervals_2 = [integer_intervals_2, interval_2];

    interval_3 = interval_1 + random_number_2; 
    interval_3(interval_3 > 719) = interval_3(interval_3 > 719) - 720;
    integer_intervals_3 = [integer_intervals_3, interval_3];

    interval_4 = interval_1 + random_number_3; 
    interval_4(interval_4 > 719) = interval_4(interval_4 > 719) - 720;
    integer_intervals_4 = [integer_intervals_4, interval_4];

    interval_5 = interval_1 + n; 
    interval_5(interval_5 > 719) = interval_5(interval_5 > 719) - 720;
    integer_intervals_5 = [integer_intervals_5, interval_5];
   
end


% 두 행렬을 합쳐서 하나의 행렬로 만듭니다.
sum_intervals = ...
    [integer_intervals_1, ...
    integer_intervals_2, ...
    integer_intervals_3, ...
    integer_intervals_4, ...
    integer_intervals_5, ...
    ];

% 행렬을 오름차순으로 정렬하고 중복된 요소를 제거합니다.
unique_intervals = unique(sort(sum_intervals));


% 두 행렬의 길이 차이 계산
difference = length(sum_intervals) - length(unique_intervals) ;

 % 결과를 작업 공간에 저장
    result_4 = [result_4; n, difference];
end

A_4 = length(unique_intervals) / 720 * 100;

% result 변수에서 2열이 최소인 행을 찾습니다.
rows_with_min_4 = result_4(result_4(:, 2) == min(result_4(:, 2)), :);

k_4 = rows_with_min_4(:,1);


% 랜덤 인덱스 선택
random_index_4 = randi([1, length(k_4)]);

% 선택된 숫자
random_number_4 = k_4(random_index_4);

%% Random 5

for n = 0:719
integer_intervals_1 = [];
integer_intervals_2 = [];
integer_intervals_3 = [];
integer_intervals_4 = [];
integer_intervals_5 = [];
integer_intervals_6 = [];


for i = 1:length(ff_0)
    interval_1 = start_points(i):end_points(i);
    interval_1(interval_1 > 719) = interval_1(interval_1 > 719) - 720; 
    integer_intervals_1 = [integer_intervals_1, interval_1];

    interval_2 = interval_1 + random_number_1; 
    interval_2(interval_2 > 719) = interval_2(interval_2 > 719) - 720;
    integer_intervals_2 = [integer_intervals_2, interval_2];

    interval_3 = interval_1 + random_number_2; 
    interval_3(interval_3 > 719) = interval_3(interval_3 > 719) - 720;
    integer_intervals_3 = [integer_intervals_3, interval_3];

    interval_4 = interval_1 + random_number_3; 
    interval_4(interval_4 > 719) = interval_4(interval_4 > 719) - 720;
    integer_intervals_4 = [integer_intervals_4, interval_4];

    interval_5 = interval_1 + random_number_4; 
    interval_5(interval_5 > 719) = interval_5(interval_5 > 719) - 720;
    integer_intervals_5 = [integer_intervals_5, interval_5];

    interval_6 = interval_1 + n; 
    interval_6(interval_6 > 719) = interval_6(interval_6 > 719) - 720;
    integer_intervals_6 = [integer_intervals_6, interval_6];
   
end


% 두 행렬을 합쳐서 하나의 행렬로 만듭니다.
sum_intervals = ...
    [integer_intervals_1, ...
    integer_intervals_2, ...
    integer_intervals_3, ...
    integer_intervals_4, ...
    integer_intervals_5, ...
    integer_intervals_6, ...
    ];

% 행렬을 오름차순으로 정렬하고 중복된 요소를 제거합니다.
unique_intervals = unique(sort(sum_intervals));


% 두 행렬의 길이 차이 계산
difference = length(sum_intervals) - length(unique_intervals) ;

 % 결과를 작업 공간에 저장
    result_5 = [result_5; n, difference];
end

A_5 = length(unique_intervals) / 720 * 100;

% result 변수에서 2열이 최소인 행을 찾습니다.
rows_with_min_5 = result_5(result_5(:, 2) == min(result_5(:, 2)), :);

k_5 = rows_with_min_5(:,1);


% 랜덤 인덱스 선택
random_index_5 = randi([1, length(k_5)]);

% 선택된 숫자
random_number_5 = k_5(random_index_5);

%% Random 6

for n = 0:719
integer_intervals_1 = [];
integer_intervals_2 = [];
integer_intervals_3 = [];
integer_intervals_4 = [];
integer_intervals_5 = [];
integer_intervals_6 = [];
integer_intervals_7 = [];


for i = 1:length(ff_0)
    interval_1 = start_points(i):end_points(i);
    interval_1(interval_1 > 719) = interval_1(interval_1 > 719) - 720; 
    integer_intervals_1 = [integer_intervals_1, interval_1];

    interval_2 = interval_1 + random_number_1; 
    interval_2(interval_2 > 719) = interval_2(interval_2 > 719) - 720;
    integer_intervals_2 = [integer_intervals_2, interval_2];

    interval_3 = interval_1 + random_number_2; 
    interval_3(interval_3 > 719) = interval_3(interval_3 > 719) - 720;
    integer_intervals_3 = [integer_intervals_3, interval_3];

    interval_4 = interval_1 + random_number_3; 
    interval_4(interval_4 > 719) = interval_4(interval_4 > 719) - 720;
    integer_intervals_4 = [integer_intervals_4, interval_4];

    interval_5 = interval_1 + random_number_4; 
    interval_5(interval_5 > 719) = interval_5(interval_5 > 719) - 720;
    integer_intervals_5 = [integer_intervals_5, interval_5];

    interval_6 = interval_1 + random_number_5; 
    interval_6(interval_6 > 719) = interval_6(interval_6 > 719) - 720;
    integer_intervals_6 = [integer_intervals_6, interval_6];

    interval_7 = interval_1 + n; 
    interval_7(interval_7 > 719) = interval_7(interval_7 > 719) - 720;
    integer_intervals_7 = [integer_intervals_7, interval_7];
   
end


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
    result_6 = [result_6; n, difference];
end

A_6 = length(unique_intervals) / 720 * 100;

% result 변수에서 2열이 최소인 행을 찾습니다.
rows_with_min_6 = result_6(result_6(:, 2) == min(result_6(:, 2)), :);

k_6 = rows_with_min_6(:,1);


% 랜덤 인덱스 선택
random_index_6 = randi([1, length(k_6)]);

% 선택된 숫자
random_number_6 = k_6(random_index_6);

%% Random 7

for n = 0:719
integer_intervals_1 = [];
integer_intervals_2 = [];
integer_intervals_3 = [];
integer_intervals_4 = [];
integer_intervals_5 = [];
integer_intervals_6 = [];
integer_intervals_7 = [];
integer_intervals_8 = [];


for i = 1:length(ff_0)
    interval_1 = start_points(i):end_points(i);
    interval_1(interval_1 > 719) = interval_1(interval_1 > 719) - 720; 
    integer_intervals_1 = [integer_intervals_1, interval_1];

    interval_2 = interval_1 + random_number_1; 
    interval_2(interval_2 > 719) = interval_2(interval_2 > 719) - 720;
    integer_intervals_2 = [integer_intervals_2, interval_2];

    interval_3 = interval_1 + random_number_2; 
    interval_3(interval_3 > 719) = interval_3(interval_3 > 719) - 720;
    integer_intervals_3 = [integer_intervals_3, interval_3];

    interval_4 = interval_1 + random_number_3; 
    interval_4(interval_4 > 719) = interval_4(interval_4 > 719) - 720;
    integer_intervals_4 = [integer_intervals_4, interval_4];

    interval_5 = interval_1 + random_number_4; 
    interval_5(interval_5 > 719) = interval_5(interval_5 > 719) - 720;
    integer_intervals_5 = [integer_intervals_5, interval_5];

    interval_6 = interval_1 + random_number_5; 
    interval_6(interval_6 > 719) = interval_6(interval_6 > 719) - 720;
    integer_intervals_6 = [integer_intervals_6, interval_6];

    interval_7 = interval_1 +  random_number_6; 
    interval_7(interval_7 > 719) = interval_7(interval_7 > 719) - 720;
    integer_intervals_7 = [integer_intervals_7, interval_7];

     interval_8 = interval_1 + n; 
    interval_8(interval_8 > 719) = interval_8(interval_8 > 719) - 720;
    integer_intervals_8 = [integer_intervals_8, interval_8];
   
end


% 두 행렬을 합쳐서 하나의 행렬로 만듭니다.
sum_intervals = ...
    [integer_intervals_1, ...
    integer_intervals_2, ...
    integer_intervals_3, ...
    integer_intervals_4, ...
    integer_intervals_5, ...
    integer_intervals_6, ...
    integer_intervals_7, ...
    integer_intervals_8, ...
    ];

% 행렬을 오름차순으로 정렬하고 중복된 요소를 제거합니다.
unique_intervals = unique(sort(sum_intervals));


% 두 행렬의 길이 차이 계산
difference = length(sum_intervals) - length(unique_intervals) ;

 % 결과를 작업 공간에 저장
    result_7 = [result_7; n, difference];
end

A_7 = length(unique_intervals) / 720 * 100;

% result 변수에서 2열이 최소인 행을 찾습니다.
rows_with_min_7 = result_7(result_7(:, 2) == min(result_7(:, 2)), :);

k_7 = rows_with_min_7(:,1);


% 랜덤 인덱스 선택
random_index_7 = randi([1, length(k_7)]);

% 선택된 숫자
random_number_7 = k_7(random_index_7);

%% Random 8

for n = 0:719
integer_intervals_1 = [];
integer_intervals_2 = [];
integer_intervals_3 = [];
integer_intervals_4 = [];
integer_intervals_5 = [];
integer_intervals_6 = [];
integer_intervals_7 = [];
integer_intervals_8 = [];
integer_intervals_9 = [];


for i = 1:length(ff_0)
    interval_1 = start_points(i):end_points(i);
    interval_1(interval_1 > 719) = interval_1(interval_1 > 719) - 720; 
    integer_intervals_1 = [integer_intervals_1, interval_1];

    interval_2 = interval_1 + random_number_1; 
    interval_2(interval_2 > 719) = interval_2(interval_2 > 719) - 720;
    integer_intervals_2 = [integer_intervals_2, interval_2];

    interval_3 = interval_1 + random_number_2; 
    interval_3(interval_3 > 719) = interval_3(interval_3 > 719) - 720;
    integer_intervals_3 = [integer_intervals_3, interval_3];

    interval_4 = interval_1 + random_number_3; 
    interval_4(interval_4 > 719) = interval_4(interval_4 > 719) - 720;
    integer_intervals_4 = [integer_intervals_4, interval_4];

    interval_5 = interval_1 + random_number_4; 
    interval_5(interval_5 > 719) = interval_5(interval_5 > 719) - 720;
    integer_intervals_5 = [integer_intervals_5, interval_5];

    interval_6 = interval_1 + random_number_5; 
    interval_6(interval_6 > 719) = interval_6(interval_6 > 719) - 720;
    integer_intervals_6 = [integer_intervals_6, interval_6];

    interval_7 = interval_1 +  random_number_6; 
    interval_7(interval_7 > 719) = interval_7(interval_7 > 719) - 720;
    integer_intervals_7 = [integer_intervals_7, interval_7];

    interval_8 = interval_1 + random_number_7; 
    interval_8(interval_8 > 719) = interval_8(interval_8 > 719) - 720;
    integer_intervals_8 = [integer_intervals_8, interval_8];
   
    interval_9 = interval_1 + n; 
    interval_9(interval_9 > 719) = interval_9(interval_9 > 719) - 720;
    integer_intervals_9 = [integer_intervals_9, interval_9];

end


% 두 행렬을 합쳐서 하나의 행렬로 만듭니다.
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
    ];

% 행렬을 오름차순으로 정렬하고 중복된 요소를 제거합니다.
unique_intervals = unique(sort(sum_intervals));


% 두 행렬의 길이 차이 계산
difference = length(sum_intervals) - length(unique_intervals) ;

 % 결과를 작업 공간에 저장
    result_8 = [result_8; n, difference];
end

A_8 = length(unique_intervals) / 720 * 100;

% result 변수에서 2열이 최소인 행을 찾습니다.
rows_with_min_8 = result_8(result_8(:, 2) == min(result_8(:, 2)), :);

k_8 = rows_with_min_8(:,1);


% 랜덤 인덱스 선택
random_index_8 = randi([1, length(k_8)]);

% 선택된 숫자
random_number_8 = k_8(random_index_8);

%% Random 9

for n = 0:719
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


for i = 1:length(ff_0)
    interval_1 = start_points(i):end_points(i);
    interval_1(interval_1 > 719) = interval_1(interval_1 > 719) - 720; 
    integer_intervals_1 = [integer_intervals_1, interval_1];

    interval_2 = interval_1 + random_number_1; 
    interval_2(interval_2 > 719) = interval_2(interval_2 > 719) - 720;
    integer_intervals_2 = [integer_intervals_2, interval_2];

    interval_3 = interval_1 + random_number_2; 
    interval_3(interval_3 > 719) = interval_3(interval_3 > 719) - 720;
    integer_intervals_3 = [integer_intervals_3, interval_3];

    interval_4 = interval_1 + random_number_3; 
    interval_4(interval_4 > 719) = interval_4(interval_4 > 719) - 720;
    integer_intervals_4 = [integer_intervals_4, interval_4];

    interval_5 = interval_1 + random_number_4; 
    interval_5(interval_5 > 719) = interval_5(interval_5 > 719) - 720;
    integer_intervals_5 = [integer_intervals_5, interval_5];

    interval_6 = interval_1 + random_number_5; 
    interval_6(interval_6 > 719) = interval_6(interval_6 > 719) - 720;
    integer_intervals_6 = [integer_intervals_6, interval_6];

    interval_7 = interval_1 +  random_number_6; 
    interval_7(interval_7 > 719) = interval_7(interval_7 > 719) - 720;
    integer_intervals_7 = [integer_intervals_7, interval_7];

    interval_8 = interval_1 + random_number_7; 
    interval_8(interval_8 > 719) = interval_8(interval_8 > 719) - 720;
    integer_intervals_8 = [integer_intervals_8, interval_8];
   
    interval_9 = interval_1 + random_number_8; 
    interval_9(interval_9 > 719) = interval_9(interval_9 > 719) - 720;
    integer_intervals_9 = [integer_intervals_9, interval_9];

    interval_10 = interval_1 + n; 
    interval_10(interval_10 > 719) = interval_10(interval_10 > 719) - 720;
    integer_intervals_10 = [integer_intervals_10, interval_10];

end


% 두 행렬을 합쳐서 하나의 행렬로 만듭니다.
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
    ];

% 행렬을 오름차순으로 정렬하고 중복된 요소를 제거합니다.
unique_intervals = unique(sort(sum_intervals));


% 두 행렬의 길이 차이 계산
difference = length(sum_intervals) - length(unique_intervals) ;

 % 결과를 작업 공간에 저장
    result_9 = [result_9; n, difference];
end

A_9 = length(unique_intervals) / 720 * 100;

% result 변수에서 2열이 최소인 행을 찾습니다.
rows_with_min_9 = result_9(result_9(:, 2) == min(result_9(:, 2)), :);

k_9 = rows_with_min_9(:,1);


% 랜덤 인덱스 선택
random_index_9 = randi([1, length(k_9)]);

% 선택된 숫자
random_number_9 = k_9(random_index_9);

%% Random 10

for n = 0:719
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


for i = 1:length(ff_0)
    interval_1 = start_points(i):end_points(i);
    interval_1(interval_1 > 719) = interval_1(interval_1 > 719) - 720; 
    integer_intervals_1 = [integer_intervals_1, interval_1];

    interval_2 = interval_1 + random_number_1; 
    interval_2(interval_2 > 719) = interval_2(interval_2 > 719) - 720;
    integer_intervals_2 = [integer_intervals_2, interval_2];

    interval_3 = interval_1 + random_number_2; 
    interval_3(interval_3 > 719) = interval_3(interval_3 > 719) - 720;
    integer_intervals_3 = [integer_intervals_3, interval_3];

    interval_4 = interval_1 + random_number_3; 
    interval_4(interval_4 > 719) = interval_4(interval_4 > 719) - 720;
    integer_intervals_4 = [integer_intervals_4, interval_4];

    interval_5 = interval_1 + random_number_4; 
    interval_5(interval_5 > 719) = interval_5(interval_5 > 719) - 720;
    integer_intervals_5 = [integer_intervals_5, interval_5];

    interval_6 = interval_1 + random_number_5; 
    interval_6(interval_6 > 719) = interval_6(interval_6 > 719) - 720;
    integer_intervals_6 = [integer_intervals_6, interval_6];

    interval_7 = interval_1 +  random_number_6; 
    interval_7(interval_7 > 719) = interval_7(interval_7 > 719) - 720;
    integer_intervals_7 = [integer_intervals_7, interval_7];

    interval_8 = interval_1 + random_number_7; 
    interval_8(interval_8 > 719) = interval_8(interval_8 > 719) - 720;
    integer_intervals_8 = [integer_intervals_8, interval_8];
   
    interval_9 = interval_1 + random_number_8; 
    interval_9(interval_9 > 719) = interval_9(interval_9 > 719) - 720;
    integer_intervals_9 = [integer_intervals_9, interval_9];

    interval_10 = interval_1 + random_number_9; 
    interval_10(interval_10 > 719) = interval_10(interval_10 > 719) - 720;
    integer_intervals_10 = [integer_intervals_10, interval_10];

    interval_11 = interval_1 + n; 
    interval_11(interval_11 > 719) = interval_11(interval_11 > 719) - 720;
    integer_intervals_11 = [integer_intervals_11, interval_11];

end


% 두 행렬을 합쳐서 하나의 행렬로 만듭니다.
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
    ];

% 행렬을 오름차순으로 정렬하고 중복된 요소를 제거합니다.
unique_intervals = unique(sort(sum_intervals));


% 두 행렬의 길이 차이 계산
difference = length(sum_intervals) - length(unique_intervals) ;

 % 결과를 작업 공간에 저장
    result_10 = [result_10; n, difference];
end

A_10 = length(unique_intervals) / 720 * 100;

% result 변수에서 2열이 최소인 행을 찾습니다.
rows_with_min_10 = result_10(result_10(:, 2) == min(result_10(:, 2)), :);

k_10 = rows_with_min_10(:,1);


% 랜덤 인덱스 선택
random_index_10 = randi([1, length(k_10)]);

% 선택된 숫자
random_number_10 = k_10(random_index_10);

%% Random 11

for n = 0:719
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


for i = 1:length(ff_0)
    interval_1 = start_points(i):end_points(i);
    interval_1(interval_1 > 719) = interval_1(interval_1 > 719) - 720; 
    integer_intervals_1 = [integer_intervals_1, interval_1];

    interval_2 = interval_1 + random_number_1; 
    interval_2(interval_2 > 719) = interval_2(interval_2 > 719) - 720;
    integer_intervals_2 = [integer_intervals_2, interval_2];

    interval_3 = interval_1 + random_number_2; 
    interval_3(interval_3 > 719) = interval_3(interval_3 > 719) - 720;
    integer_intervals_3 = [integer_intervals_3, interval_3];

    interval_4 = interval_1 + random_number_3; 
    interval_4(interval_4 > 719) = interval_4(interval_4 > 719) - 720;
    integer_intervals_4 = [integer_intervals_4, interval_4];

    interval_5 = interval_1 + random_number_4; 
    interval_5(interval_5 > 719) = interval_5(interval_5 > 719) - 720;
    integer_intervals_5 = [integer_intervals_5, interval_5];

    interval_6 = interval_1 + random_number_5; 
    interval_6(interval_6 > 719) = interval_6(interval_6 > 719) - 720;
    integer_intervals_6 = [integer_intervals_6, interval_6];

    interval_7 = interval_1 +  random_number_6; 
    interval_7(interval_7 > 719) = interval_7(interval_7 > 719) - 720;
    integer_intervals_7 = [integer_intervals_7, interval_7];

    interval_8 = interval_1 + random_number_7; 
    interval_8(interval_8 > 719) = interval_8(interval_8 > 719) - 720;
    integer_intervals_8 = [integer_intervals_8, interval_8];
   
    interval_9 = interval_1 + random_number_8; 
    interval_9(interval_9 > 719) = interval_9(interval_9 > 719) - 720;
    integer_intervals_9 = [integer_intervals_9, interval_9];

    interval_10 = interval_1 + random_number_9; 
    interval_10(interval_10 > 719) = interval_10(interval_10 > 719) - 720;
    integer_intervals_10 = [integer_intervals_10, interval_10];

    interval_11 = interval_1 + random_number_10; 
    interval_11(interval_11 > 719) = interval_11(interval_11 > 719) - 720;
    integer_intervals_11 = [integer_intervals_11, interval_11];

    interval_12 = interval_1 + n; 
    interval_12(interval_12 > 719) = interval_12(interval_12 > 719) - 720;
    integer_intervals_12 = [integer_intervals_12, interval_12];

end


% 두 행렬을 합쳐서 하나의 행렬로 만듭니다.
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
    ];

% 행렬을 오름차순으로 정렬하고 중복된 요소를 제거합니다.
unique_intervals = unique(sort(sum_intervals));


% 두 행렬의 길이 차이 계산
difference = length(sum_intervals) - length(unique_intervals) ;

 % 결과를 작업 공간에 저장
    result_11 = [result_11; n, difference];
end

A_11 = length(unique_intervals) / 720 * 100;

% result 변수에서 2열이 최소인 행을 찾습니다.
rows_with_min_11 = result_11(result_11(:, 2) == min(result_11(:, 2)), :);

k_11 = rows_with_min_11(:,1);


% 랜덤 인덱스 선택
random_index_11 = randi([1, length(k_11)]);

% 선택된 숫자
random_number_11 = k_11(random_index_11);

%% Random 12

for n = 0:719
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


for i = 1:length(ff_0)
    interval_1 = start_points(i):end_points(i);
    interval_1(interval_1 > 719) = interval_1(interval_1 > 719) - 720; 
    integer_intervals_1 = [integer_intervals_1, interval_1];

    interval_2 = interval_1 + random_number_1; 
    interval_2(interval_2 > 719) = interval_2(interval_2 > 719) - 720;
    integer_intervals_2 = [integer_intervals_2, interval_2];

    interval_3 = interval_1 + random_number_2; 
    interval_3(interval_3 > 719) = interval_3(interval_3 > 719) - 720;
    integer_intervals_3 = [integer_intervals_3, interval_3];

    interval_4 = interval_1 + random_number_3; 
    interval_4(interval_4 > 719) = interval_4(interval_4 > 719) - 720;
    integer_intervals_4 = [integer_intervals_4, interval_4];

    interval_5 = interval_1 + random_number_4; 
    interval_5(interval_5 > 719) = interval_5(interval_5 > 719) - 720;
    integer_intervals_5 = [integer_intervals_5, interval_5];

    interval_6 = interval_1 + random_number_5; 
    interval_6(interval_6 > 719) = interval_6(interval_6 > 719) - 720;
    integer_intervals_6 = [integer_intervals_6, interval_6];

    interval_7 = interval_1 +  random_number_6; 
    interval_7(interval_7 > 719) = interval_7(interval_7 > 719) - 720;
    integer_intervals_7 = [integer_intervals_7, interval_7];

    interval_8 = interval_1 + random_number_7; 
    interval_8(interval_8 > 719) = interval_8(interval_8 > 719) - 720;
    integer_intervals_8 = [integer_intervals_8, interval_8];
   
    interval_9 = interval_1 + random_number_8; 
    interval_9(interval_9 > 719) = interval_9(interval_9 > 719) - 720;
    integer_intervals_9 = [integer_intervals_9, interval_9];

    interval_10 = interval_1 + random_number_9; 
    interval_10(interval_10 > 719) = interval_10(interval_10 > 719) - 720;
    integer_intervals_10 = [integer_intervals_10, interval_10];

    interval_11 = interval_1 + random_number_10; 
    interval_11(interval_11 > 719) = interval_11(interval_11 > 719) - 720;
    integer_intervals_11 = [integer_intervals_11, interval_11];

    interval_12 = interval_1 + random_number_11; 
    interval_12(interval_12 > 719) = interval_12(interval_12 > 719) - 720;
    integer_intervals_12 = [integer_intervals_12, interval_12];

    interval_13 = interval_1 + n; 
    interval_13(interval_13 > 719) = interval_13(interval_13 > 719) - 720;
    integer_intervals_13 = [integer_intervals_13, interval_13];

end


% 두 행렬을 합쳐서 하나의 행렬로 만듭니다.
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
    ];

% 행렬을 오름차순으로 정렬하고 중복된 요소를 제거합니다.
unique_intervals = unique(sort(sum_intervals));


% 두 행렬의 길이 차이 계산
difference = length(sum_intervals) - length(unique_intervals) ;

 % 결과를 작업 공간에 저장
    result_12 = [result_12; n, difference];
end

A_12 = length(unique_intervals) / 720 * 100;

% result 변수에서 2열이 최소인 행을 찾습니다.
rows_with_min_12 = result_12(result_12(:, 2) == min(result_12(:, 2)), :);

k_12 = rows_with_min_12(:,1);


% 랜덤 인덱스 선택
random_index_12= randi([1, length(k_12)]);

% 선택된 숫자
random_number_12 = k_12(random_index_12);

%% Random 13

for n = 0:719
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


for i = 1:length(ff_0)
    interval_1 = start_points(i):end_points(i);
    interval_1(interval_1 > 719) = interval_1(interval_1 > 719) - 720; 
    integer_intervals_1 = [integer_intervals_1, interval_1];

    interval_2 = interval_1 + random_number_1; 
    interval_2(interval_2 > 719) = interval_2(interval_2 > 719) - 720;
    integer_intervals_2 = [integer_intervals_2, interval_2];

    interval_3 = interval_1 + random_number_2; 
    interval_3(interval_3 > 719) = interval_3(interval_3 > 719) - 720;
    integer_intervals_3 = [integer_intervals_3, interval_3];

    interval_4 = interval_1 + random_number_3; 
    interval_4(interval_4 > 719) = interval_4(interval_4 > 719) - 720;
    integer_intervals_4 = [integer_intervals_4, interval_4];

    interval_5 = interval_1 + random_number_4; 
    interval_5(interval_5 > 719) = interval_5(interval_5 > 719) - 720;
    integer_intervals_5 = [integer_intervals_5, interval_5];

    interval_6 = interval_1 + random_number_5; 
    interval_6(interval_6 > 719) = interval_6(interval_6 > 719) - 720;
    integer_intervals_6 = [integer_intervals_6, interval_6];

    interval_7 = interval_1 +  random_number_6; 
    interval_7(interval_7 > 719) = interval_7(interval_7 > 719) - 720;
    integer_intervals_7 = [integer_intervals_7, interval_7];

    interval_8 = interval_1 + random_number_7; 
    interval_8(interval_8 > 719) = interval_8(interval_8 > 719) - 720;
    integer_intervals_8 = [integer_intervals_8, interval_8];
   
    interval_9 = interval_1 + random_number_8; 
    interval_9(interval_9 > 719) = interval_9(interval_9 > 719) - 720;
    integer_intervals_9 = [integer_intervals_9, interval_9];

    interval_10 = interval_1 + random_number_9; 
    interval_10(interval_10 > 719) = interval_10(interval_10 > 719) - 720;
    integer_intervals_10 = [integer_intervals_10, interval_10];

    interval_11 = interval_1 + random_number_10; 
    interval_11(interval_11 > 719) = interval_11(interval_11 > 719) - 720;
    integer_intervals_11 = [integer_intervals_11, interval_11];

    interval_12 = interval_1 + random_number_11; 
    interval_12(interval_12 > 719) = interval_12(interval_12 > 719) - 720;
    integer_intervals_12 = [integer_intervals_12, interval_12];

    interval_13 = interval_1 + random_number_12; 
    interval_13(interval_13 > 719) = interval_13(interval_13 > 719) - 720;
    integer_intervals_13 = [integer_intervals_13, interval_13];

    interval_14 = interval_1 + n; 
    interval_14(interval_14 > 719) = interval_14(interval_14 > 719) - 720;
    integer_intervals_14 = [integer_intervals_14, interval_14];

end


% 두 행렬을 합쳐서 하나의 행렬로 만듭니다.
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
    ];

% 행렬을 오름차순으로 정렬하고 중복된 요소를 제거합니다.
unique_intervals = unique(sort(sum_intervals));


% 두 행렬의 길이 차이 계산
difference = length(sum_intervals) - length(unique_intervals) ;

 % 결과를 작업 공간에 저장
    result_13 = [result_13; n, difference];
end

A_13 = length(unique_intervals) / 720 * 100;

% result 변수에서 2열이 최소인 행을 찾습니다.
rows_with_min_13 = result_13(result_13(:, 2) == min(result_13(:, 2)), :);

k_13 = rows_with_min_13(:,1);


% 랜덤 인덱스 선택
random_index_13= randi([1, length(k_13)]);

% 선택된 숫자
random_number_13 = k_13(random_index_13);

%% Random 14

for n = 0:719
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


for i = 1:length(ff_0)
    interval_1 = start_points(i):end_points(i);
    interval_1(interval_1 > 719) = interval_1(interval_1 > 719) - 720; 
    integer_intervals_1 = [integer_intervals_1, interval_1];

    interval_2 = interval_1 + random_number_1; 
    interval_2(interval_2 > 719) = interval_2(interval_2 > 719) - 720;
    integer_intervals_2 = [integer_intervals_2, interval_2];

    interval_3 = interval_1 + random_number_2; 
    interval_3(interval_3 > 719) = interval_3(interval_3 > 719) - 720;
    integer_intervals_3 = [integer_intervals_3, interval_3];

    interval_4 = interval_1 + random_number_3; 
    interval_4(interval_4 > 719) = interval_4(interval_4 > 719) - 720;
    integer_intervals_4 = [integer_intervals_4, interval_4];

    interval_5 = interval_1 + random_number_4; 
    interval_5(interval_5 > 719) = interval_5(interval_5 > 719) - 720;
    integer_intervals_5 = [integer_intervals_5, interval_5];

    interval_6 = interval_1 + random_number_5; 
    interval_6(interval_6 > 719) = interval_6(interval_6 > 719) - 720;
    integer_intervals_6 = [integer_intervals_6, interval_6];

    interval_7 = interval_1 +  random_number_6; 
    interval_7(interval_7 > 719) = interval_7(interval_7 > 719) - 720;
    integer_intervals_7 = [integer_intervals_7, interval_7];

    interval_8 = interval_1 + random_number_7; 
    interval_8(interval_8 > 719) = interval_8(interval_8 > 719) - 720;
    integer_intervals_8 = [integer_intervals_8, interval_8];
   
    interval_9 = interval_1 + random_number_8; 
    interval_9(interval_9 > 719) = interval_9(interval_9 > 719) - 720;
    integer_intervals_9 = [integer_intervals_9, interval_9];

    interval_10 = interval_1 + random_number_9; 
    interval_10(interval_10 > 719) = interval_10(interval_10 > 719) - 720;
    integer_intervals_10 = [integer_intervals_10, interval_10];

    interval_11 = interval_1 + random_number_10; 
    interval_11(interval_11 > 719) = interval_11(interval_11 > 719) - 720;
    integer_intervals_11 = [integer_intervals_11, interval_11];

    interval_12 = interval_1 + random_number_11; 
    interval_12(interval_12 > 719) = interval_12(interval_12 > 719) - 720;
    integer_intervals_12 = [integer_intervals_12, interval_12];

    interval_13 = interval_1 + random_number_12; 
    interval_13(interval_13 > 719) = interval_13(interval_13 > 719) - 720;
    integer_intervals_13 = [integer_intervals_13, interval_13];

    interval_14 = interval_1 + random_number_13; 
    interval_14(interval_14 > 719) = interval_14(interval_14 > 719) - 720;
    integer_intervals_14 = [integer_intervals_14, interval_14];

    interval_15 = interval_1 + n; 
    interval_15(interval_15 > 719) = interval_15(interval_15 > 719) - 720;
    integer_intervals_15 = [integer_intervals_15, interval_15];


end


% 두 행렬을 합쳐서 하나의 행렬로 만듭니다.
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
    ];

% 행렬을 오름차순으로 정렬하고 중복된 요소를 제거합니다.
unique_intervals = unique(sort(sum_intervals));


% 두 행렬의 길이 차이 계산
difference = length(sum_intervals) - length(unique_intervals) ;

 % 결과를 작업 공간에 저장
    result_14 = [result_14; n, difference];
end

A_14 = length(unique_intervals) / 720 * 100;

% result 변수에서 2열이 최소인 행을 찾습니다.
rows_with_min_14 = result_14(result_14(:, 2) == min(result_14(:, 2)), :);

k_14 = rows_with_min_14(:,1);


% 랜덤 인덱스 선택
random_index_14= randi([1, length(k_14)]);

% 선택된 숫자
random_number_14 = k_14(random_index_14);

%% Random 15

for n = 0:719
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


for i = 1:length(ff_0)
    interval_1 = start_points(i):end_points(i);
    interval_1(interval_1 > 719) = interval_1(interval_1 > 719) - 720; 
    integer_intervals_1 = [integer_intervals_1, interval_1];

    interval_2 = interval_1 + random_number_1; 
    interval_2(interval_2 > 719) = interval_2(interval_2 > 719) - 720;
    integer_intervals_2 = [integer_intervals_2, interval_2];

    interval_3 = interval_1 + random_number_2; 
    interval_3(interval_3 > 719) = interval_3(interval_3 > 719) - 720;
    integer_intervals_3 = [integer_intervals_3, interval_3];

    interval_4 = interval_1 + random_number_3; 
    interval_4(interval_4 > 719) = interval_4(interval_4 > 719) - 720;
    integer_intervals_4 = [integer_intervals_4, interval_4];

    interval_5 = interval_1 + random_number_4; 
    interval_5(interval_5 > 719) = interval_5(interval_5 > 719) - 720;
    integer_intervals_5 = [integer_intervals_5, interval_5];

    interval_6 = interval_1 + random_number_5; 
    interval_6(interval_6 > 719) = interval_6(interval_6 > 719) - 720;
    integer_intervals_6 = [integer_intervals_6, interval_6];

    interval_7 = interval_1 +  random_number_6; 
    interval_7(interval_7 > 719) = interval_7(interval_7 > 719) - 720;
    integer_intervals_7 = [integer_intervals_7, interval_7];

    interval_8 = interval_1 + random_number_7; 
    interval_8(interval_8 > 719) = interval_8(interval_8 > 719) - 720;
    integer_intervals_8 = [integer_intervals_8, interval_8];
   
    interval_9 = interval_1 + random_number_8; 
    interval_9(interval_9 > 719) = interval_9(interval_9 > 719) - 720;
    integer_intervals_9 = [integer_intervals_9, interval_9];

    interval_10 = interval_1 + random_number_9; 
    interval_10(interval_10 > 719) = interval_10(interval_10 > 719) - 720;
    integer_intervals_10 = [integer_intervals_10, interval_10];

    interval_11 = interval_1 + random_number_10; 
    interval_11(interval_11 > 719) = interval_11(interval_11 > 719) - 720;
    integer_intervals_11 = [integer_intervals_11, interval_11];

    interval_12 = interval_1 + random_number_11; 
    interval_12(interval_12 > 719) = interval_12(interval_12 > 719) - 720;
    integer_intervals_12 = [integer_intervals_12, interval_12];

    interval_13 = interval_1 + random_number_12; 
    interval_13(interval_13 > 719) = interval_13(interval_13 > 719) - 720;
    integer_intervals_13 = [integer_intervals_13, interval_13];

    interval_14 = interval_1 + random_number_13; 
    interval_14(interval_14 > 719) = interval_14(interval_14 > 719) - 720;
    integer_intervals_14 = [integer_intervals_14, interval_14];

    interval_15 = interval_1 + random_number_14; 
    interval_15(interval_15 > 719) = interval_15(interval_15 > 719) - 720;
    integer_intervals_15 = [integer_intervals_15, interval_15];

    interval_16 = interval_1 + n; 
    interval_16(interval_16 > 719) = interval_16(interval_16 > 719) - 720;
    integer_intervals_16 = [integer_intervals_16, interval_16];


end


% 두 행렬을 합쳐서 하나의 행렬로 만듭니다.
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
    ];

% 행렬을 오름차순으로 정렬하고 중복된 요소를 제거합니다.
unique_intervals = unique(sort(sum_intervals));


% 두 행렬의 길이 차이 계산
difference = length(sum_intervals) - length(unique_intervals) ;

 % 결과를 작업 공간에 저장
    result_15 = [result_15; n, difference];
end

A_15 = length(unique_intervals) / 720 * 100;

% result 변수에서 2열이 최소인 행을 찾습니다.
rows_with_min_15 = result_15(result_15(:, 2) == min(result_15(:, 2)), :);

k_15 = rows_with_min_15(:,1);


% 랜덤 인덱스 선택
random_index_15= randi([1, length(k_15)]);

% 선택된 숫자
random_number_15 = k_15(random_index_15);

%% Random 16

for n = 0:719
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


for i = 1:length(ff_0)
    interval_1 = start_points(i):end_points(i);
    interval_1(interval_1 > 719) = interval_1(interval_1 > 719) - 720; 
    integer_intervals_1 = [integer_intervals_1, interval_1];

    interval_2 = interval_1 + random_number_1; 
    interval_2(interval_2 > 719) = interval_2(interval_2 > 719) - 720;
    integer_intervals_2 = [integer_intervals_2, interval_2];

    interval_3 = interval_1 + random_number_2; 
    interval_3(interval_3 > 719) = interval_3(interval_3 > 719) - 720;
    integer_intervals_3 = [integer_intervals_3, interval_3];

    interval_4 = interval_1 + random_number_3; 
    interval_4(interval_4 > 719) = interval_4(interval_4 > 719) - 720;
    integer_intervals_4 = [integer_intervals_4, interval_4];

    interval_5 = interval_1 + random_number_4; 
    interval_5(interval_5 > 719) = interval_5(interval_5 > 719) - 720;
    integer_intervals_5 = [integer_intervals_5, interval_5];

    interval_6 = interval_1 + random_number_5; 
    interval_6(interval_6 > 719) = interval_6(interval_6 > 719) - 720;
    integer_intervals_6 = [integer_intervals_6, interval_6];

    interval_7 = interval_1 +  random_number_6; 
    interval_7(interval_7 > 719) = interval_7(interval_7 > 719) - 720;
    integer_intervals_7 = [integer_intervals_7, interval_7];

    interval_8 = interval_1 + random_number_7; 
    interval_8(interval_8 > 719) = interval_8(interval_8 > 719) - 720;
    integer_intervals_8 = [integer_intervals_8, interval_8];
   
    interval_9 = interval_1 + random_number_8; 
    interval_9(interval_9 > 719) = interval_9(interval_9 > 719) - 720;
    integer_intervals_9 = [integer_intervals_9, interval_9];

    interval_10 = interval_1 + random_number_9; 
    interval_10(interval_10 > 719) = interval_10(interval_10 > 719) - 720;
    integer_intervals_10 = [integer_intervals_10, interval_10];

    interval_11 = interval_1 + random_number_10; 
    interval_11(interval_11 > 719) = interval_11(interval_11 > 719) - 720;
    integer_intervals_11 = [integer_intervals_11, interval_11];

    interval_12 = interval_1 + random_number_11; 
    interval_12(interval_12 > 719) = interval_12(interval_12 > 719) - 720;
    integer_intervals_12 = [integer_intervals_12, interval_12];

    interval_13 = interval_1 + random_number_12; 
    interval_13(interval_13 > 719) = interval_13(interval_13 > 719) - 720;
    integer_intervals_13 = [integer_intervals_13, interval_13];

    interval_14 = interval_1 + random_number_13; 
    interval_14(interval_14 > 719) = interval_14(interval_14 > 719) - 720;
    integer_intervals_14 = [integer_intervals_14, interval_14];

    interval_15 = interval_1 + random_number_14; 
    interval_15(interval_15 > 719) = interval_15(interval_15 > 719) - 720;
    integer_intervals_15 = [integer_intervals_15, interval_15];

    interval_16 = interval_1 + random_number_15; 
    interval_16(interval_16 > 719) = interval_16(interval_16 > 719) - 720;
    integer_intervals_16 = [integer_intervals_16, interval_16];

    interval_17 = interval_1 + n; 
    interval_17(interval_17 > 719) = interval_17(interval_17 > 719) - 720;
    integer_intervals_17 = [integer_intervals_17, interval_17];


end


% 두 행렬을 합쳐서 하나의 행렬로 만듭니다.
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
    ];

% 행렬을 오름차순으로 정렬하고 중복된 요소를 제거합니다.
unique_intervals = unique(sort(sum_intervals));


% 두 행렬의 길이 차이 계산
difference = length(sum_intervals) - length(unique_intervals) ;

 % 결과를 작업 공간에 저장
    result_16 = [result_16; n, difference];
end

A_16 = length(unique_intervals) / 720 * 100;

% result 변수에서 2열이 최소인 행을 찾습니다.
rows_with_min_16 = result_16(result_16(:, 2) == min(result_16(:, 2)), :);

k_16 = rows_with_min_16(:,1);


% 랜덤 인덱스 선택
random_index_16= randi([1, length(k_16)]);

% 선택된 숫자
random_number_16 = k_16(random_index_16);

%% Random 17

for n = 0:719
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


for i = 1:length(ff_0)
    interval_1 = start_points(i):end_points(i);
    interval_1(interval_1 > 719) = interval_1(interval_1 > 719) - 720; 
    integer_intervals_1 = [integer_intervals_1, interval_1];

    interval_2 = interval_1 + random_number_1; 
    interval_2(interval_2 > 719) = interval_2(interval_2 > 719) - 720;
    integer_intervals_2 = [integer_intervals_2, interval_2];

    interval_3 = interval_1 + random_number_2; 
    interval_3(interval_3 > 719) = interval_3(interval_3 > 719) - 720;
    integer_intervals_3 = [integer_intervals_3, interval_3];

    interval_4 = interval_1 + random_number_3; 
    interval_4(interval_4 > 719) = interval_4(interval_4 > 719) - 720;
    integer_intervals_4 = [integer_intervals_4, interval_4];

    interval_5 = interval_1 + random_number_4; 
    interval_5(interval_5 > 719) = interval_5(interval_5 > 719) - 720;
    integer_intervals_5 = [integer_intervals_5, interval_5];

    interval_6 = interval_1 + random_number_5; 
    interval_6(interval_6 > 719) = interval_6(interval_6 > 719) - 720;
    integer_intervals_6 = [integer_intervals_6, interval_6];

    interval_7 = interval_1 +  random_number_6; 
    interval_7(interval_7 > 719) = interval_7(interval_7 > 719) - 720;
    integer_intervals_7 = [integer_intervals_7, interval_7];

    interval_8 = interval_1 + random_number_7; 
    interval_8(interval_8 > 719) = interval_8(interval_8 > 719) - 720;
    integer_intervals_8 = [integer_intervals_8, interval_8];
   
    interval_9 = interval_1 + random_number_8; 
    interval_9(interval_9 > 719) = interval_9(interval_9 > 719) - 720;
    integer_intervals_9 = [integer_intervals_9, interval_9];

    interval_10 = interval_1 + random_number_9; 
    interval_10(interval_10 > 719) = interval_10(interval_10 > 719) - 720;
    integer_intervals_10 = [integer_intervals_10, interval_10];

    interval_11 = interval_1 + random_number_10; 
    interval_11(interval_11 > 719) = interval_11(interval_11 > 719) - 720;
    integer_intervals_11 = [integer_intervals_11, interval_11];

    interval_12 = interval_1 + random_number_11; 
    interval_12(interval_12 > 719) = interval_12(interval_12 > 719) - 720;
    integer_intervals_12 = [integer_intervals_12, interval_12];

    interval_13 = interval_1 + random_number_12; 
    interval_13(interval_13 > 719) = interval_13(interval_13 > 719) - 720;
    integer_intervals_13 = [integer_intervals_13, interval_13];

    interval_14 = interval_1 + random_number_13; 
    interval_14(interval_14 > 719) = interval_14(interval_14 > 719) - 720;
    integer_intervals_14 = [integer_intervals_14, interval_14];

    interval_15 = interval_1 + random_number_14; 
    interval_15(interval_15 > 719) = interval_15(interval_15 > 719) - 720;
    integer_intervals_15 = [integer_intervals_15, interval_15];

    interval_16 = interval_1 + random_number_15; 
    interval_16(interval_16 > 719) = interval_16(interval_16 > 719) - 720;
    integer_intervals_16 = [integer_intervals_16, interval_16];

    interval_17 = interval_1 + random_number_16; 
    interval_17(interval_17 > 719) = interval_17(interval_17 > 719) - 720;
    integer_intervals_17 = [integer_intervals_17, interval_17];

    interval_18 = interval_1 + n; 
    interval_18(interval_18 > 719) = interval_18(interval_18 > 719) - 720;
    integer_intervals_18 = [integer_intervals_18, interval_18];


end


% 두 행렬을 합쳐서 하나의 행렬로 만듭니다.
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
    ];

% 행렬을 오름차순으로 정렬하고 중복된 요소를 제거합니다.
unique_intervals = unique(sort(sum_intervals));


% 두 행렬의 길이 차이 계산
difference = length(sum_intervals) - length(unique_intervals) ;

 % 결과를 작업 공간에 저장
    result_17 = [result_17; n, difference];
end

A_17 = length(unique_intervals) / 720 * 100;

% result 변수에서 2열이 최소인 행을 찾습니다.
rows_with_min_17 = result_17(result_17(:, 2) == min(result_17(:, 2)), :);

k_17 = rows_with_min_17(:,1);


% 랜덤 인덱스 선택
random_index_17= randi([1, length(k_17)]);

% 선택된 숫자
random_number_17 = k_17(random_index_17);

%% Random 18

for n = 0:719
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


for i = 1:length(ff_0)
    interval_1 = start_points(i):end_points(i);
    interval_1(interval_1 > 719) = interval_1(interval_1 > 719) - 720; 
    integer_intervals_1 = [integer_intervals_1, interval_1];

    interval_2 = interval_1 + random_number_1; 
    interval_2(interval_2 > 719) = interval_2(interval_2 > 719) - 720;
    integer_intervals_2 = [integer_intervals_2, interval_2];

    interval_3 = interval_1 + random_number_2; 
    interval_3(interval_3 > 719) = interval_3(interval_3 > 719) - 720;
    integer_intervals_3 = [integer_intervals_3, interval_3];

    interval_4 = interval_1 + random_number_3; 
    interval_4(interval_4 > 719) = interval_4(interval_4 > 719) - 720;
    integer_intervals_4 = [integer_intervals_4, interval_4];

    interval_5 = interval_1 + random_number_4; 
    interval_5(interval_5 > 719) = interval_5(interval_5 > 719) - 720;
    integer_intervals_5 = [integer_intervals_5, interval_5];

    interval_6 = interval_1 + random_number_5; 
    interval_6(interval_6 > 719) = interval_6(interval_6 > 719) - 720;
    integer_intervals_6 = [integer_intervals_6, interval_6];

    interval_7 = interval_1 +  random_number_6; 
    interval_7(interval_7 > 719) = interval_7(interval_7 > 719) - 720;
    integer_intervals_7 = [integer_intervals_7, interval_7];

    interval_8 = interval_1 + random_number_7; 
    interval_8(interval_8 > 719) = interval_8(interval_8 > 719) - 720;
    integer_intervals_8 = [integer_intervals_8, interval_8];
   
    interval_9 = interval_1 + random_number_8; 
    interval_9(interval_9 > 719) = interval_9(interval_9 > 719) - 720;
    integer_intervals_9 = [integer_intervals_9, interval_9];

    interval_10 = interval_1 + random_number_9; 
    interval_10(interval_10 > 719) = interval_10(interval_10 > 719) - 720;
    integer_intervals_10 = [integer_intervals_10, interval_10];

    interval_11 = interval_1 + random_number_10; 
    interval_11(interval_11 > 719) = interval_11(interval_11 > 719) - 720;
    integer_intervals_11 = [integer_intervals_11, interval_11];

    interval_12 = interval_1 + random_number_11; 
    interval_12(interval_12 > 719) = interval_12(interval_12 > 719) - 720;
    integer_intervals_12 = [integer_intervals_12, interval_12];

    interval_13 = interval_1 + random_number_12; 
    interval_13(interval_13 > 719) = interval_13(interval_13 > 719) - 720;
    integer_intervals_13 = [integer_intervals_13, interval_13];

    interval_14 = interval_1 + random_number_13; 
    interval_14(interval_14 > 719) = interval_14(interval_14 > 719) - 720;
    integer_intervals_14 = [integer_intervals_14, interval_14];

    interval_15 = interval_1 + random_number_14; 
    interval_15(interval_15 > 719) = interval_15(interval_15 > 719) - 720;
    integer_intervals_15 = [integer_intervals_15, interval_15];

    interval_16 = interval_1 + random_number_15; 
    interval_16(interval_16 > 719) = interval_16(interval_16 > 719) - 720;
    integer_intervals_16 = [integer_intervals_16, interval_16];

    interval_17 = interval_1 + random_number_16; 
    interval_17(interval_17 > 719) = interval_17(interval_17 > 719) - 720;
    integer_intervals_17 = [integer_intervals_17, interval_17];

    interval_18 = interval_1 + random_number_17; 
    interval_18(interval_18 > 719) = interval_18(interval_18 > 719) - 720;
    integer_intervals_18 = [integer_intervals_18, interval_18];

    interval_19 = interval_1 + n; 
    interval_19(interval_19 > 719) = interval_19(interval_19 > 719) - 720;
    integer_intervals_19 = [integer_intervals_19, interval_19];


end


% 두 행렬을 합쳐서 하나의 행렬로 만듭니다.
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
    ];

% 행렬을 오름차순으로 정렬하고 중복된 요소를 제거합니다.
unique_intervals = unique(sort(sum_intervals));


% 두 행렬의 길이 차이 계산
difference = length(sum_intervals) - length(unique_intervals) ;

 % 결과를 작업 공간에 저장
    result_18 = [result_18; n, difference];
end

A_18 = length(unique_intervals) / 720 * 100

% result 변수에서 2열이 최소인 행을 찾습니다.
rows_with_min_18 = result_18(result_18(:, 2) == min(result_18(:, 2)), :);

k_18 = rows_with_min_18(:,1);

 


random_N(1) = (random_number_1);
random_N(2) = (random_number_2);
random_N(3) = (random_number_3);
random_N(4) = (random_number_4);
random_N(5) = (random_number_5);
random_N(6) = (random_number_6);
random_N(7) = (random_number_7);
random_N(8) = (random_number_8);
random_N(9) = (random_number_9);
random_N(10) = (random_number_10);
random_N(11) = (random_number_11);
random_N(12) = (random_number_12);
random_N(13) = (random_number_13);
random_N(14) = (random_number_14);
random_N(15) = (random_number_15);
random_N(16) = (random_number_16);
random_N(17) = (random_number_17)

  if A_18 >= 99


        break;
    end
 trials = trials + 1; 
end

% 시뮬레이션 종료 시간 기록
elapsed_time = toc(start_time);

% 결과 출력
fprintf('시행 횟수: %d\n', trials);
fprintf('걸린 시간(초): %.2f\n', elapsed_time);

