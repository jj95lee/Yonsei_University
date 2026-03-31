clc
clear all

A_7 = 0;
trials = 0;
start_time = tic;
L = 717;
shift = 0;

while true

% 주어진 행렬 ff_0
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

% spot1
% 1422.34375	3882.988281
% 12744.29688	15186.19141
% 23907.98828	26012.96875
% 35150.78125	35906.23047
% 55517.20703	56159.45313
% 65349.64844	67430.07813
% 76151.81641	78589.90234

% spot 2
% 3639.394531	6092.890625
% 42206.21094	44808.82813
% 53656.60156	55843.61328
% 78515	80951.67969



% ff_0 -> n 으로 변환
ff_0 = floor(ff_00/(86024/L)) - shift;
neg1 = ff_0;

length_ap = ff_0(1,2)-ff_0(1,1)+ff_0(2,2)-ff_0(2,1)+ff_0(3,2)-ff_0(3,1)+...
+ff_0(4,2)-ff_0(4,1)+4;



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




%% Random 1

for n = 0:716
integer_intervals_1 = [];
integer_intervals_2 = [];


for i = 1:length(ff_0)
    interval_1 = start_points(i) :end_points(i) ;
    interval_1(interval_1 > 716) = interval_1(interval_1 > 716) - 717; 
    integer_intervals_1 = [integer_intervals_1, interval_1];

    interval_2 = interval_1 + n; 
    interval_2(interval_2 > 716) = interval_2(interval_2 > 716) - 717;
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

A_1 = length(unique_intervals) / 717 * 100;

% result 변수에서 2열이 최소인 행을 찾습니다.
rows_with_min_1 = result_1(result_1(:, 2) == min(result_1(:, 2)), :);

k_1 = rows_with_min_1(:,1);


% 랜덤 인덱스 선택
random_index_1 = randi([1, length(k_1)]);

% 선택된 숫자
random_number_1 = k_1(random_index_1);



%% Random 2

for n = 0:716
integer_intervals_1 = [];
integer_intervals_2 = [];
integer_intervals_3 = [];


for i = 1:length(ff_0)
    interval_1 = start_points(i) :end_points(i);
    interval_1(interval_1 > 716) = interval_1(interval_1 > 716) - 717; 
    integer_intervals_1 = [integer_intervals_1, interval_1];

    interval_2 = interval_1 + random_number_1; 
    interval_2(interval_2 > 716) = interval_2(interval_2 > 716) - 717;
    integer_intervals_2 = [integer_intervals_2, interval_2];

    interval_3 = interval_1 + n; 
    interval_3(interval_3 > 716) = interval_3(interval_3 > 716) - 717;
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

A_2 = length(unique_intervals) / 717 * 100;

% result 변수에서 2열이 최소인 행을 찾습니다.
rows_with_min_2 = result_2(result_2(:, 2) == min(result_2(:, 2)), :);

k_2 = rows_with_min_2(:,1);


% 랜덤 인덱스 선택
random_index_2 = randi([1, length(k_2)]);

% 선택된 숫자
random_number_2 = k_2(random_index_2);


%% Random 3

for n = 0:716
integer_intervals_1 = [];
integer_intervals_2 = [];
integer_intervals_3 = [];
integer_intervals_4 = [];


for i = 1:length(ff_0)
    interval_1 = start_points(i):end_points(i);
    interval_1(interval_1 > 716) = interval_1(interval_1 > 716) - 717; 
    integer_intervals_1 = [integer_intervals_1, interval_1];

    interval_2 = interval_1 + random_number_1; 
    interval_2(interval_2 > 716) = interval_2(interval_2 > 716) - 717;
    integer_intervals_2 = [integer_intervals_2, interval_2];

    interval_3 = interval_1 + random_number_2; 
    interval_3(interval_3 > 716) = interval_3(interval_3 > 716) - 717;
    integer_intervals_3 = [integer_intervals_3, interval_3];

    interval_4 = interval_1 + n; 
    interval_4(interval_4 > 716) = interval_4(interval_4 > 716) - 717;
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

A_3 = length(unique_intervals) / 717 * 100;

% result 변수에서 2열이 최소인 행을 찾습니다.
rows_with_min_3 = result_3(result_3(:, 2) == min(result_3(:, 2)), :);

k_3 = rows_with_min_3(:,1);


% 랜덤 인덱스 선택
random_index_3 = randi([1, length(k_3)]);

% 선택된 숫자
random_number_3 = k_3(random_index_3);

%% Random 4

for n = 0:716
integer_intervals_1 = [];
integer_intervals_2 = [];
integer_intervals_3 = [];
integer_intervals_4 = [];
integer_intervals_5 = [];


for i = 1:length(ff_0)
    interval_1 = start_points(i):end_points(i);
    interval_1(interval_1 > 716) = interval_1(interval_1 > 716) - 717; 
    integer_intervals_1 = [integer_intervals_1, interval_1];

    interval_2 = interval_1 + random_number_1; 
    interval_2(interval_2 > 716) = interval_2(interval_2 > 716) - 717;
    integer_intervals_2 = [integer_intervals_2, interval_2];

    interval_3 = interval_1 + random_number_2; 
    interval_3(interval_3 > 716) = interval_3(interval_3 > 716) - 717;
    integer_intervals_3 = [integer_intervals_3, interval_3];

    interval_4 = interval_1 + random_number_3; 
    interval_4(interval_4 > 716) = interval_4(interval_4 > 716) - 717;
    integer_intervals_4 = [integer_intervals_4, interval_4];

    interval_5 = interval_1 + n; 
    interval_5(interval_5 > 716) = interval_5(interval_5 > 716) - 717;
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

A_4 = length(unique_intervals) / 717 * 100;

% result 변수에서 2열이 최소인 행을 찾습니다.
rows_with_min_4 = result_4(result_4(:, 2) == min(result_4(:, 2)), :);

k_4 = rows_with_min_4(:,1);


% 랜덤 인덱스 선택
random_index_4 = randi([1, length(k_4)]);

% 선택된 숫자
random_number_4 = k_4(random_index_4);

%% Random 5

for n = 0:716
integer_intervals_1 = [];
integer_intervals_2 = [];
integer_intervals_3 = [];
integer_intervals_4 = [];
integer_intervals_5 = [];
integer_intervals_6 = [];


for i = 1:length(ff_0)
    interval_1 = start_points(i):end_points(i);
    interval_1(interval_1 > 716) = interval_1(interval_1 > 716) - 717; 
    integer_intervals_1 = [integer_intervals_1, interval_1];

    interval_2 = interval_1 + random_number_1; 
    interval_2(interval_2 > 716) = interval_2(interval_2 > 716) - 717;
    integer_intervals_2 = [integer_intervals_2, interval_2];

    interval_3 = interval_1 + random_number_2; 
    interval_3(interval_3 > 716) = interval_3(interval_3 > 716) - 717;
    integer_intervals_3 = [integer_intervals_3, interval_3];

    interval_4 = interval_1 + random_number_3; 
    interval_4(interval_4 > 716) = interval_4(interval_4 > 716) - 717;
    integer_intervals_4 = [integer_intervals_4, interval_4];

    interval_5 = interval_1 + random_number_4; 
    interval_5(interval_5 > 716) = interval_5(interval_5 > 716) - 717;
    integer_intervals_5 = [integer_intervals_5, interval_5];

    interval_6 = interval_1 + n; 
    interval_6(interval_6 > 716) = interval_6(interval_6 > 716) - 717;
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

A_5 = length(unique_intervals) / 717 * 100;

% result 변수에서 2열이 최소인 행을 찾습니다.
rows_with_min_5 = result_5(result_5(:, 2) == min(result_5(:, 2)), :);

k_5 = rows_with_min_5(:,1);


% 랜덤 인덱스 선택
random_index_5 = randi([1, length(k_5)]);

% 선택된 숫자
random_number_5 = k_5(random_index_5);

%% Random 6

for n = 0:716
integer_intervals_1 = [];
integer_intervals_2 = [];
integer_intervals_3 = [];
integer_intervals_4 = [];
integer_intervals_5 = [];
integer_intervals_6 = [];
integer_intervals_7 = [];


for i = 1:length(ff_0)
    interval_1 = start_points(i):end_points(i);
    interval_1(interval_1 > 716) = interval_1(interval_1 > 716) - 717; 
    integer_intervals_1 = [integer_intervals_1, interval_1];

    interval_2 = interval_1 + random_number_1; 
    interval_2(interval_2 > 716) = interval_2(interval_2 > 716) - 717;
    integer_intervals_2 = [integer_intervals_2, interval_2];

    interval_3 = interval_1 + random_number_2; 
    interval_3(interval_3 > 716) = interval_3(interval_3 > 716) - 717;
    integer_intervals_3 = [integer_intervals_3, interval_3];

    interval_4 = interval_1 + random_number_3; 
    interval_4(interval_4 > 716) = interval_4(interval_4 > 716) - 717;
    integer_intervals_4 = [integer_intervals_4, interval_4];

    interval_5 = interval_1 + random_number_4; 
    interval_5(interval_5 > 716) = interval_5(interval_5 > 716) - 717;
    integer_intervals_5 = [integer_intervals_5, interval_5];

    interval_6 = interval_1 + random_number_5; 
    interval_6(interval_6 > 716) = interval_6(interval_6 > 716) - 717;
    integer_intervals_6 = [integer_intervals_6, interval_6];

    interval_7 = interval_1 + n; 
    interval_7(interval_7 > 716) = interval_7(interval_7 > 716) - 717;
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

A_6 = length(unique_intervals) / 717 * 100;

% result 변수에서 2열이 최소인 행을 찾습니다.
rows_with_min_6 = result_6(result_6(:, 2) == min(result_6(:, 2)), :);

k_6 = rows_with_min_6(:,1);


% 랜덤 인덱스 선택
random_index_6 = randi([1, length(k_6)]);

% 선택된 숫자
random_number_6 = k_6(random_index_6);


%% Random 7

for n = 0:716
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
    interval_1(interval_1 > 716) = interval_1(interval_1 > 716) - 717; 
    integer_intervals_1 = [integer_intervals_1, interval_1];

    interval_2 = interval_1 + random_number_1; 
    interval_2(interval_2 > 716) = interval_2(interval_2 > 716) - 717;
    integer_intervals_2 = [integer_intervals_2, interval_2];

    interval_3 = interval_1 + random_number_2; 
    interval_3(interval_3 > 716) = interval_3(interval_3 > 716) - 717;
    integer_intervals_3 = [integer_intervals_3, interval_3];

    interval_4 = interval_1 + random_number_3; 
    interval_4(interval_4 > 716) = interval_4(interval_4 > 716) - 717;
    integer_intervals_4 = [integer_intervals_4, interval_4];

    interval_5 = interval_1 + random_number_4; 
    interval_5(interval_5 > 716) = interval_5(interval_5 > 716) - 717;
    integer_intervals_5 = [integer_intervals_5, interval_5];

    interval_6 = interval_1 + random_number_5; 
    interval_6(interval_6 > 716) = interval_6(interval_6 > 716) - 717;
    integer_intervals_6 = [integer_intervals_6, interval_6];

    interval_7 = interval_1 +  random_number_6; 
    interval_7(interval_7 > 716) = interval_7(interval_7 > 716) - 717;
    integer_intervals_7 = [integer_intervals_7, interval_7];

    interval_8 = interval_1 + n; 
    interval_8(interval_8 > 716) = interval_8(interval_8 > 716) - 717;
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
    ];

% 행렬을 오름차순으로 정렬하고 중복된 요소를 제거합니다.
unique_intervals = unique(sort(sum_intervals));


% 두 행렬의 길이 차이 계산
difference = length(sum_intervals) - length(unique_intervals) ;

 % 결과를 작업 공간에 저장
    result_7 = [result_7; n, difference];
end

A_7 = length(unique_intervals) / 717 * 100;

% result 변수에서 2열이 최소인 행을 찾습니다.
rows_with_min_7 = result_7(result_7(:, 2) == min(result_7(:, 2)), :);

k_7 = rows_with_min_7(:,1);

% 랜덤 인덱스 선택
random_index_7 = randi([1, length(k_7)]);

% 선택된 숫자
random_number_7 = k_7(random_index_7);


%% Random 8

for n = 0:716
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
    interval_1(interval_1 > 716) = interval_1(interval_1 > 716) - 717; 
    integer_intervals_1 = [integer_intervals_1, interval_1];

    interval_2 = interval_1 + random_number_1; 
    interval_2(interval_2 > 716) = interval_2(interval_2 > 716) - 717;
    integer_intervals_2 = [integer_intervals_2, interval_2];

    interval_3 = interval_1 + random_number_2; 
    interval_3(interval_3 > 716) = interval_3(interval_3 > 716) - 717;
    integer_intervals_3 = [integer_intervals_3, interval_3];

    interval_4 = interval_1 + random_number_3; 
    interval_4(interval_4 > 716) = interval_4(interval_4 > 716) - 717;
    integer_intervals_4 = [integer_intervals_4, interval_4];

    interval_5 = interval_1 + random_number_4; 
    interval_5(interval_5 > 716) = interval_5(interval_5 > 716) - 717;
    integer_intervals_5 = [integer_intervals_5, interval_5];

    interval_6 = interval_1 + random_number_5; 
    interval_6(interval_6 > 716) = interval_6(interval_6 > 716) - 717;
    integer_intervals_6 = [integer_intervals_6, interval_6];

    interval_7 = interval_1 +  random_number_6; 
    interval_7(interval_7 > 716) = interval_7(interval_7 > 716) - 717;
    integer_intervals_7 = [integer_intervals_7, interval_7];

    interval_8 = interval_1 + random_number_7; 
    interval_8(interval_8 > 716) = interval_8(interval_8 > 716) - 717;
    integer_intervals_8 = [integer_intervals_8, interval_8];

    interval_9 = interval_1 + n; 
    interval_8(interval_8 > 716) = interval_8(interval_8 > 716) - 717;
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
    result_8 = [result_8; n, difference];
end

A_8 = length(unique_intervals) / 717 * 100

% result 변수에서 2열이 최소인 행을 찾습니다.
rows_with_min_8 = result_8(result_8(:, 2) == min(result_8(:, 2)), :);

k_8 = rows_with_min_8(:,1);

% 랜덤 인덱스 선택
random_index_8 = randi([1, length(k_8)]);

% 선택된 숫자
random_number_8 = k_8(random_index_8);

%% Random 9

for n = 0:716
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
    interval_1(interval_1 > 716) = interval_1(interval_1 > 716) - 717; 
    integer_intervals_1 = [integer_intervals_1, interval_1];

    interval_2 = interval_1 + random_number_1; 
    interval_2(interval_2 > 716) = interval_2(interval_2 > 716) - 717;
    integer_intervals_2 = [integer_intervals_2, interval_2];

    interval_3 = interval_1 + random_number_2; 
    interval_3(interval_3 > 716) = interval_3(interval_3 > 716) - 717;
    integer_intervals_3 = [integer_intervals_3, interval_3];

    interval_4 = interval_1 + random_number_3; 
    interval_4(interval_4 > 716) = interval_4(interval_4 > 716) - 717;
    integer_intervals_4 = [integer_intervals_4, interval_4];

    interval_5 = interval_1 + random_number_4; 
    interval_5(interval_5 > 716) = interval_5(interval_5 > 716) - 717;
    integer_intervals_5 = [integer_intervals_5, interval_5];

    interval_6 = interval_1 + random_number_5; 
    interval_6(interval_6 > 716) = interval_6(interval_6 > 716) - 717;
    integer_intervals_6 = [integer_intervals_6, interval_6];

    interval_7 = interval_1 +  random_number_6; 
    interval_7(interval_7 > 716) = interval_7(interval_7 > 716) - 717;
    integer_intervals_7 = [integer_intervals_7, interval_7];

    interval_8 = interval_1 + random_number_7; 
    interval_8(interval_8 > 716) = interval_8(interval_8 > 716) - 717;
    integer_intervals_8 = [integer_intervals_8, interval_8];

    interval_9 = interval_1 + random_number_8; 
    interval_8(interval_8 > 716) = interval_8(interval_8 > 716) - 717;
    integer_intervals_8 = [integer_intervals_8, interval_8];

    interval_10 = interval_1 + n; 
    interval_9(interval_9 > 716) = interval_9(interval_9 > 716) - 717;
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
    result_9 = [result_9; n, difference];
end

A_9 = length(unique_intervals) / 717 * 100

% result 변수에서 2열이 최소인 행을 찾습니다.
rows_with_min_9 = result_9(result_9(:, 2) == min(result_9(:, 2)), :);

k_9 = rows_with_min_9(:,1);

% 랜덤 인덱스 선택
random_index_9 = randi([1, length(k_9)]);

% 선택된 숫자
random_number_9 = k_9(random_index_9);
 
%%

random_N(1) = (random_number_1);
random_N(2) = (random_number_2);
random_N(3) = (random_number_3);
random_N(4) = (random_number_4);
random_N(5) = (random_number_5);
random_N(6) = (random_number_6);
random_N(7) = (random_number_7);
random_N(8) = (random_number_8)

  if A_9 >= 99.9


        break;
    end
 trials = trials + 1; 
end

% 시뮬레이션 종료 시간 기록
elapsed_time = toc(start_time);

% 결과 출력
fprintf('시행 횟수: %d\n', trials);
fprintf('걸린 시간(초): %.2f\n', elapsed_time);

