clear all

trials = 0;
start_time = tic;
L = 720;
shift = 0;

while true

% 주어진 행렬 ff_0
ff_00 = ...
[...
40897.26563	41339.23828
46910.625	47431.28906
52987.67578	53507.51953
59062.55859	59578.82813
65174.88281	65556.21094
...		
];

gg_00 = ...
[...
35018.67188	35214.08203
40899.375	41407.91016
46972.96875	47481.79688
53072.8125	53572.85156
59146.75781	59667.1875
65260.37109	65652.94922
...		
];


% 저궤도

% spot1 korea
% 40897.26563	41339.23828
% 46910.625	47431.28906
% 52987.67578	53507.51953
% 59062.55859	59578.82813
% 65174.88281	65556.21094


% spot 2 example
% 35018.67188	35214.08203
% 40899.375	41407.91016
% 46972.96875	47481.79688
% 53072.8125	53572.85156
% 59146.75781	59667.1875
% 65260.37109	65652.94922



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
ff_0 = floor(ff_00/(86400/L)) - shift;
neg1 = ff_0;

gg_0 = floor(gg_00/(86400/L)) - shift;
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


%% Seed Satellite Setting(0 / 507~512)

for n = 0:719
integer_intervals_1 = [];
integer_intervals_2 = [];
integer_intervals_3 = [];

   % 저궤도 seed satellite 설정
   for i = 1:length(ff_0)
       interval_1 = start_points1(i) :end_points1(i) ;
       interval_1(interval_1 > 719) = interval_1(interval_1 > 719) - 720; 
       integer_intervals_1 = [integer_intervals_1, interval_1];
   end

   % 고궤도 seed satellite 틀만 설정
   for i = 1:length(gg_0)
       interval_2 = start_points2(i) :end_points2(i) ;
       interval_2(interval_2 > 719) = interval_2(interval_2 > 719) - 720; 
       integer_intervals_2 = [integer_intervals_2, interval_2];

   % 고궤도 seed satellite 랜덤 위치 생성
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

   % 고궤도 seed satellite 틀만 설정
   for i = 1:length(gg_0)
       interval_2 = start_points2(i) :end_points2(i) ;
       interval_2(interval_2 > 719) = interval_2(interval_2 > 719) - 720; 
       integer_intervals_2 = [integer_intervals_2, interval_2];

       % 고궤도 seed satellite 위치 부여(507~512)
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


%% 저궤도 추가 (2L / 1H)

for n = 0:719
integer_intervals_1 = [];
integer_intervals_2 = [];
integer_intervals_3 = [];
integer_intervals_4 = [];


for i = 1:length(ff_0)
    interval_1 = start_points1(i) :end_points1(i) ;
    interval_1(interval_1 > 719) = interval_1(interval_1 > 719) - 720; 
    integer_intervals_1 = [integer_intervals_1, interval_1];
end


for i = 1:length(gg_0)
    interval_2 = start_points2(i) :end_points2(i) ;
    interval_2(interval_2 > 719) = interval_2(interval_2 > 719) - 720; 
    integer_intervals_2 = [integer_intervals_2, interval_2];

    interval_3 = interval_2 + random_number_1; 
    interval_3(interval_3 > 719) = interval_3(interval_3 > 719) - 720;
    integer_intervals_3 = [integer_intervals_3, interval_3];
end

for i = 1:length(ff_0)
    interval_1 = start_points1(i) :end_points1(i) ;
    interval_1(interval_1 > 719) = interval_1(interval_1 > 719) - 720; 
    integer_intervals_1 = [integer_intervals_1, interval_1];

    interval_4 = interval_1 + n; 
    interval_4(interval_4 > 719) = interval_4(interval_4 > 719) - 720;
    integer_intervals_4 = [integer_intervals_4, interval_4];
end


% 두 행렬을 합쳐서 하나의 행렬로 만듭니다.
sum_intervals = ...
    [integer_intervals_1, ...
    integer_intervals_3, ...
    integer_intervals_4, ...
    ];

% 행렬을 오름차순으로 정렬하고 중복된 요소를 제거합니다.
unique_intervals = unique(sort(sum_intervals));


% 두 행렬의 길이 차이 계산
difference = length(sum_intervals) - length(unique_intervals) ;

 % 결과를 작업 공간에 저장
    result_2 = [result_2; n, difference];
end


% result 변수에서 2열이 최소인 행을 찾습니다.
rows_with_min_2 = result_2(result_2(:, 2) == min(result_2(:, 2)), :);

k_2 = rows_with_min_2(:,1);


% 랜덤 인덱스 선택
random_index_2 = randi([1, length(k_2)]);

% 선택된 숫자
random_number_2 = k_2(random_index_2);

% new
for n = 0:719
integer_intervals_1 = [];
integer_intervals_2 = [];
integer_intervals_3 = [];
integer_intervals_4 = [];


for i = 1:length(ff_0)
    interval_1 = start_points1(i) :end_points1(i) ;
    interval_1(interval_1 > 719) = interval_1(interval_1 > 719) - 720; 
    integer_intervals_1 = [integer_intervals_1, interval_1];
end


for i = 1:length(gg_0)
    interval_2 = start_points2(i) :end_points2(i) ;
    interval_2(interval_2 > 719) = interval_2(interval_2 > 719) - 720; 
    integer_intervals_2 = [integer_intervals_2, interval_2];

    interval_3 = interval_2 + random_number_1; 
    interval_3(interval_3 > 719) = interval_3(interval_3 > 719) - 720;
    integer_intervals_3 = [integer_intervals_3, interval_3];
end

for i = 1:length(ff_0)
    interval_1 = start_points1(i) :end_points1(i) ;
    interval_1(interval_1 > 719) = interval_1(interval_1 > 719) - 720; 
    integer_intervals_1 = [integer_intervals_1, interval_1];

    interval_4 = interval_1 + random_number_2; 
    interval_4(interval_4 > 719) = interval_4(interval_4 > 719) - 720;
    integer_intervals_4 = [integer_intervals_4, interval_4];
end

% 두 행렬을 합쳐서 하나의 행렬로 만듭니다.
sum_intervals = ...
    [integer_intervals_1, ...
    integer_intervals_3, ...
    integer_intervals_4, ...
    ];

% 행렬을 오름차순으로 정렬하고 중복된 요소를 제거합니다.
unique_intervals = unique(sort(sum_intervals));

% 두 행렬의 길이 차이 계산
difference = length(sum_intervals) - length(unique_intervals) ;

 % 결과를 작업 공간에 저장
    result_2 = [result_2; n, difference];

end

A_2 = length(unique_intervals) / 720 * 100;



%% 고궤도 추가 (2L / 2H)

for n = 0:719
integer_intervals_1 = [];
integer_intervals_2 = [];
integer_intervals_3 = [];
integer_intervals_4 = [];
integer_intervals_5 = [];


for i = 1:length(ff_0)
    interval_1 = start_points1(i) :end_points1(i) ;
    interval_1(interval_1 > 719) = interval_1(interval_1 > 719) - 720; 
    integer_intervals_1 = [integer_intervals_1, interval_1];
end


for i = 1:length(gg_0)
    interval_2 = start_points2(i) :end_points2(i) ;
    interval_2(interval_2 > 719) = interval_2(interval_2 > 719) - 720; 
    integer_intervals_2 = [integer_intervals_2, interval_2];

    interval_3 = interval_2 + random_number_1; 
    interval_3(interval_3 > 719) = interval_3(interval_3 > 719) - 720;
    integer_intervals_3 = [integer_intervals_3, interval_3];
end

for i = 1:length(ff_0)
    interval_1 = start_points1(i) :end_points1(i) ;
    interval_1(interval_1 > 719) = interval_1(interval_1 > 719) - 720; 
    integer_intervals_1 = [integer_intervals_1, interval_1];

    interval_4 = interval_1 + random_number_2; 
    interval_4(interval_4 > 719) = interval_4(interval_4 > 719) - 720;
    integer_intervals_4 = [integer_intervals_4, interval_4];
end

for i = 1:length(gg_0)
    interval_2 = start_points2(i) :end_points2(i) ;
    interval_2(interval_2 > 719) = interval_2(interval_2 > 719) - 720; 
    integer_intervals_2 = [integer_intervals_2, interval_2];

    interval_5 = interval_2 + n; 
    interval_5(interval_5 > 719) = interval_5(interval_5 > 719) - 720;
    integer_intervals_5 = [integer_intervals_5, interval_5];
end


% 두 행렬을 합쳐서 하나의 행렬로 만듭니다.
sum_intervals = ...
    [integer_intervals_1, ...
    integer_intervals_3, ...
    integer_intervals_4, ...
    integer_intervals_5, ...
    ];

% 행렬을 오름차순으로 정렬하고 중복된 요소를 제거합니다.
unique_intervals = unique(sort(sum_intervals));


% 두 행렬의 길이 차이 계산
difference = length(sum_intervals) - length(unique_intervals) ;

 % 결과를 작업 공간에 저장
    result_3 = [result_3; n, difference];
end


% result 변수에서 2열이 최소인 행을 찾습니다.
rows_with_min_3 = result_3(result_3(:, 2) == min(result_3(:, 2)), :);

k_3 = rows_with_min_3(:,1);


% 랜덤 인덱스 선택
random_index_3 = randi([1, length(k_3)]);

% 선택된 숫자
random_number_3 = k_3(random_index_3);

% new
for n = 0:719
integer_intervals_1 = [];
integer_intervals_2 = [];
integer_intervals_3 = [];
integer_intervals_4 = [];
integer_intervals_5 = [];


for i = 1:length(ff_0)
    interval_1 = start_points1(i) :end_points1(i) ;
    interval_1(interval_1 > 719) = interval_1(interval_1 > 719) - 720; 
    integer_intervals_1 = [integer_intervals_1, interval_1];
end


for i = 1:length(gg_0)
    interval_2 = start_points2(i) :end_points2(i) ;
    interval_2(interval_2 > 719) = interval_2(interval_2 > 719) - 720; 
    integer_intervals_2 = [integer_intervals_2, interval_2];

    interval_3 = interval_2 + random_number_1; 
    interval_3(interval_3 > 719) = interval_3(interval_3 > 719) - 720;
    integer_intervals_3 = [integer_intervals_3, interval_3];
end

for i = 1:length(ff_0)
    interval_1 = start_points1(i) :end_points1(i) ;
    interval_1(interval_1 > 719) = interval_1(interval_1 > 719) - 720; 
    integer_intervals_1 = [integer_intervals_1, interval_1];

    interval_4 = interval_1 + random_number_2; 
    interval_4(interval_4 > 719) = interval_4(interval_4 > 719) - 720;
    integer_intervals_4 = [integer_intervals_4, interval_4];
end

for i = 1:length(gg_0)
    interval_2 = start_points2(i) :end_points2(i) ;
    interval_2(interval_2 > 719) = interval_2(interval_2 > 719) - 720; 
    integer_intervals_2 = [integer_intervals_2, interval_2];

    interval_5 = interval_2 + random_number_3; 
    interval_5(interval_5 > 719) = interval_5(interval_5 > 719) - 720;
    integer_intervals_5 = [integer_intervals_5, interval_5];
end


% 두 행렬을 합쳐서 하나의 행렬로 만듭니다.
sum_intervals = ...
    [integer_intervals_1, ...
    integer_intervals_3, ...
    integer_intervals_4, ...
    integer_intervals_5, ...
    ];

% 행렬을 오름차순으로 정렬하고 중복된 요소를 제거합니다.
unique_intervals = unique(sort(sum_intervals));


% 두 행렬의 길이 차이 계산
difference = length(sum_intervals) - length(unique_intervals) ;

 % 결과를 작업 공간에 저장
    result_3 = [result_3; n, difference];
end

A_3 = length(unique_intervals) / 720 * 100;


%% 저궤도 추가 (3L / 2H)

for n = 0:719
integer_intervals_1 = [];
integer_intervals_2 = [];
integer_intervals_3 = [];
integer_intervals_4 = [];
integer_intervals_5 = [];
integer_intervals_6 = [];


for i = 1:length(ff_0)
    interval_1 = start_points1(i) :end_points1(i) ;
    interval_1(interval_1 > 719) = interval_1(interval_1 > 719) - 720; 
    integer_intervals_1 = [integer_intervals_1, interval_1];
end


for i = 1:length(gg_0)
    interval_2 = start_points2(i) :end_points2(i) ;
    interval_2(interval_2 > 719) = interval_2(interval_2 > 719) - 720; 
    integer_intervals_2 = [integer_intervals_2, interval_2];

    interval_3 = interval_2 + random_number_1; 
    interval_3(interval_3 > 719) = interval_3(interval_3 > 719) - 720;
    integer_intervals_3 = [integer_intervals_3, interval_3];
end

for i = 1:length(ff_0)
    interval_1 = start_points1(i) :end_points1(i) ;
    interval_1(interval_1 > 719) = interval_1(interval_1 > 719) - 720; 
    integer_intervals_1 = [integer_intervals_1, interval_1];

    interval_4 = interval_1 + random_number_2; 
    interval_4(interval_4 > 719) = interval_4(interval_4 > 719) - 720;
    integer_intervals_4 = [integer_intervals_4, interval_4];
end

for i = 1:length(gg_0)
    interval_2 = start_points2(i) :end_points2(i) ;
    interval_2(interval_2 > 719) = interval_2(interval_2 > 719) - 720; 
    integer_intervals_2 = [integer_intervals_2, interval_2];

    interval_5 = interval_2 + random_number_3; 
    interval_5(interval_5 > 719) = interval_5(interval_5 > 719) - 720;
    integer_intervals_5 = [integer_intervals_5, interval_5];
end

for i = 1:length(ff_0)
    interval_1 = start_points1(i) :end_points1(i) ;
    interval_1(interval_1 > 719) = interval_1(interval_1 > 719) - 720; 
    integer_intervals_1 = [integer_intervals_1, interval_1];

    interval_6 = interval_1 + n; 
    interval_6(interval_6 > 719) = interval_6(interval_6 > 719) - 720;
    integer_intervals_6 = [integer_intervals_6, interval_6];
end


% 두 행렬을 합쳐서 하나의 행렬로 만듭니다.
sum_intervals = ...
    [integer_intervals_1, ...
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
    result_4 = [result_4; n, difference];
end

% result 변수에서 2열이 최소인 행을 찾습니다.
rows_with_min_4 = result_4(result_4(:, 2) == min(result_4(:, 2)), :);

k_4 = rows_with_min_4(:,1);


% 랜덤 인덱스 선택
random_index_4 = randi([1, length(k_4)]);

% 선택된 숫자
random_number_4 = k_4(random_index_4);

for n = 0:719
integer_intervals_1 = [];
integer_intervals_2 = [];
integer_intervals_3 = [];
integer_intervals_4 = [];
integer_intervals_5 = [];
integer_intervals_6 = [];


for i = 1:length(ff_0)
    interval_1 = start_points1(i) :end_points1(i) ;
    interval_1(interval_1 > 719) = interval_1(interval_1 > 719) - 720; 
    integer_intervals_1 = [integer_intervals_1, interval_1];
end


for i = 1:length(gg_0)
    interval_2 = start_points2(i) :end_points2(i) ;
    interval_2(interval_2 > 719) = interval_2(interval_2 > 719) - 720; 
    integer_intervals_2 = [integer_intervals_2, interval_2];

    interval_3 = interval_2 + random_number_1; 
    interval_3(interval_3 > 719) = interval_3(interval_3 > 719) - 720;
    integer_intervals_3 = [integer_intervals_3, interval_3];
end

for i = 1:length(ff_0)
    interval_1 = start_points1(i) :end_points1(i) ;
    interval_1(interval_1 > 719) = interval_1(interval_1 > 719) - 720; 
    integer_intervals_1 = [integer_intervals_1, interval_1];

    interval_4 = interval_1 + random_number_2; 
    interval_4(interval_4 > 719) = interval_4(interval_4 > 719) - 720;
    integer_intervals_4 = [integer_intervals_4, interval_4];
end

for i = 1:length(gg_0)
    interval_2 = start_points2(i) :end_points2(i) ;
    interval_2(interval_2 > 719) = interval_2(interval_2 > 719) - 720; 
    integer_intervals_2 = [integer_intervals_2, interval_2];

    interval_5 = interval_2 + random_number_3; 
    interval_5(interval_5 > 719) = interval_5(interval_5 > 719) - 720;
    integer_intervals_5 = [integer_intervals_5, interval_5];
end

for i = 1:length(ff_0)
    interval_1 = start_points1(i) :end_points1(i) ;
    interval_1(interval_1 > 719) = interval_1(interval_1 > 719) - 720; 
    integer_intervals_1 = [integer_intervals_1, interval_1];

    interval_6 = interval_1 + random_number_4; 
    interval_6(interval_6 > 719) = interval_6(interval_6 > 719) - 720;
    integer_intervals_6 = [integer_intervals_6, interval_6];
end

% 두 행렬을 합쳐서 하나의 행렬로 만듭니다.
sum_intervals = ...
    [integer_intervals_1, ...
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
    result_4 = [result_4; n, difference];
end

A_4 = length(unique_intervals) / 720 * 100;


%% 고궤도 추가 (3L / 3H)

for n = 0:719
integer_intervals_1 = [];
integer_intervals_2 = [];
integer_intervals_3 = [];
integer_intervals_4 = [];
integer_intervals_5 = [];
integer_intervals_6 = [];
integer_intervals_7 = [];


for i = 1:length(ff_0)
    interval_1 = start_points1(i) :end_points1(i) ;
    interval_1(interval_1 > 719) = interval_1(interval_1 > 719) - 720; 
    integer_intervals_1 = [integer_intervals_1, interval_1];
end


for i = 1:length(gg_0)
    interval_2 = start_points2(i) :end_points2(i) ;
    interval_2(interval_2 > 719) = interval_2(interval_2 > 719) - 720; 
    integer_intervals_2 = [integer_intervals_2, interval_2];

    interval_3 = interval_2 + random_number_1; 
    interval_3(interval_3 > 719) = interval_3(interval_3 > 719) - 720;
    integer_intervals_3 = [integer_intervals_3, interval_3];
end

for i = 1:length(ff_0)
    interval_1 = start_points1(i) :end_points1(i) ;
    interval_1(interval_1 > 719) = interval_1(interval_1 > 719) - 720; 
    integer_intervals_1 = [integer_intervals_1, interval_1];

    interval_4 = interval_1 + random_number_2; 
    interval_4(interval_4 > 719) = interval_4(interval_4 > 719) - 720;
    integer_intervals_4 = [integer_intervals_4, interval_4];
end

for i = 1:length(gg_0)
    interval_2 = start_points2(i) :end_points2(i) ;
    interval_2(interval_2 > 719) = interval_2(interval_2 > 719) - 720; 
    integer_intervals_2 = [integer_intervals_2, interval_2];

    interval_5 = interval_2 + random_number_3; 
    interval_5(interval_5 > 719) = interval_5(interval_5 > 719) - 720;
    integer_intervals_5 = [integer_intervals_5, interval_5];
end

for i = 1:length(ff_0)
    interval_1 = start_points1(i) :end_points1(i) ;
    interval_1(interval_1 > 719) = interval_1(interval_1 > 719) - 720; 
    integer_intervals_1 = [integer_intervals_1, interval_1];

    interval_6 = interval_1 + random_number_4; 
    interval_6(interval_6 > 719) = interval_6(interval_6 > 719) - 720;
    integer_intervals_6 = [integer_intervals_6, interval_6];
end

for i = 1:length(gg_0)
    interval_2 = start_points2(i) :end_points2(i) ;
    interval_2(interval_2 > 719) = interval_2(interval_2 > 719) - 720; 
    integer_intervals_2 = [integer_intervals_2, interval_2];

    interval_7 = interval_2 + n; 
    interval_7(interval_7 > 719) = interval_7(interval_7 > 719) - 720;
    integer_intervals_7 = [integer_intervals_7, interval_7];
end


% 두 행렬을 합쳐서 하나의 행렬로 만듭니다.
sum_intervals = ...
    [integer_intervals_1, ...
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
    result_5 = [result_5; n, difference];
end


% result 변수에서 2열이 최소인 행을 찾습니다.
rows_with_min_5 = result_5(result_5(:, 2) == min(result_5(:, 2)), :);

k_5 = rows_with_min_5(:,1);


% 랜덤 인덱스 선택
random_index_5 = randi([1, length(k_5)]);

% 선택된 숫자
random_number_5 = k_5(random_index_5);

for n = 0:719
integer_intervals_1 = [];
integer_intervals_2 = [];
integer_intervals_3 = [];
integer_intervals_4 = [];
integer_intervals_5 = [];
integer_intervals_6 = [];
integer_intervals_7 = [];


for i = 1:length(ff_0)
    interval_1 = start_points1(i) :end_points1(i) ;
    interval_1(interval_1 > 719) = interval_1(interval_1 > 719) - 720; 
    integer_intervals_1 = [integer_intervals_1, interval_1];
end


for i = 1:length(gg_0)
    interval_2 = start_points2(i) :end_points2(i) ;
    interval_2(interval_2 > 719) = interval_2(interval_2 > 719) - 720; 
    integer_intervals_2 = [integer_intervals_2, interval_2];

    interval_3 = interval_2 + random_number_1; 
    interval_3(interval_3 > 719) = interval_3(interval_3 > 719) - 720;
    integer_intervals_3 = [integer_intervals_3, interval_3];
end

for i = 1:length(ff_0)
    interval_1 = start_points1(i) :end_points1(i) ;
    interval_1(interval_1 > 719) = interval_1(interval_1 > 719) - 720; 
    integer_intervals_1 = [integer_intervals_1, interval_1];

    interval_4 = interval_1 + random_number_2; 
    interval_4(interval_4 > 719) = interval_4(interval_4 > 719) - 720;
    integer_intervals_4 = [integer_intervals_4, interval_4];
end

for i = 1:length(gg_0)
    interval_2 = start_points2(i) :end_points2(i) ;
    interval_2(interval_2 > 719) = interval_2(interval_2 > 719) - 720; 
    integer_intervals_2 = [integer_intervals_2, interval_2];

    interval_5 = interval_2 + random_number_3; 
    interval_5(interval_5 > 719) = interval_5(interval_5 > 719) - 720;
    integer_intervals_5 = [integer_intervals_5, interval_5];
end

for i = 1:length(ff_0)
    interval_1 = start_points1(i) :end_points1(i) ;
    interval_1(interval_1 > 719) = interval_1(interval_1 > 719) - 720; 
    integer_intervals_1 = [integer_intervals_1, interval_1];

    interval_6 = interval_1 + random_number_4; 
    interval_6(interval_6 > 719) = interval_6(interval_6 > 719) - 720;
    integer_intervals_6 = [integer_intervals_6, interval_6];
end

for i = 1:length(gg_0)
    interval_2 = start_points2(i) :end_points2(i) ;
    interval_2(interval_2 > 719) = interval_2(interval_2 > 719) - 720; 
    integer_intervals_2 = [integer_intervals_2, interval_2];

    interval_7 = interval_2 + random_number_5; 
    interval_7(interval_7 > 719) = interval_7(interval_7 > 719) - 720;
    integer_intervals_7 = [integer_intervals_7, interval_7];
end


% 두 행렬을 합쳐서 하나의 행렬로 만듭니다.
sum_intervals = ...
    [integer_intervals_1, ...
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
    result_5 = [result_5; n, difference];
end

A_5 = length(unique_intervals) / 720 * 100;



%% 저궤도 추가 (4L / 3H)

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
    interval_1 = start_points1(i) :end_points1(i) ;
    interval_1(interval_1 > 719) = interval_1(interval_1 > 719) - 720; 
    integer_intervals_1 = [integer_intervals_1, interval_1];
end


for i = 1:length(gg_0)
    interval_2 = start_points2(i) :end_points2(i) ;
    interval_2(interval_2 > 719) = interval_2(interval_2 > 719) - 720; 
    integer_intervals_2 = [integer_intervals_2, interval_2];

    interval_3 = interval_2 + random_number_1; 
    interval_3(interval_3 > 719) = interval_3(interval_3 > 719) - 720;
    integer_intervals_3 = [integer_intervals_3, interval_3];
end

for i = 1:length(ff_0)
    interval_1 = start_points1(i) :end_points1(i) ;
    interval_1(interval_1 > 719) = interval_1(interval_1 > 719) - 720; 
    integer_intervals_1 = [integer_intervals_1, interval_1];

    interval_4 = interval_1 + random_number_2; 
    interval_4(interval_4 > 719) = interval_4(interval_4 > 719) - 720;
    integer_intervals_4 = [integer_intervals_4, interval_4];
end

for i = 1:length(gg_0)
    interval_2 = start_points2(i) :end_points2(i) ;
    interval_2(interval_2 > 719) = interval_2(interval_2 > 719) - 720; 
    integer_intervals_2 = [integer_intervals_2, interval_2];

    interval_5 = interval_2 + random_number_3; 
    interval_5(interval_5 > 719) = interval_5(interval_5 > 719) - 720;
    integer_intervals_5 = [integer_intervals_5, interval_5];
end

for i = 1:length(ff_0)
    interval_1 = start_points1(i) :end_points1(i) ;
    interval_1(interval_1 > 719) = interval_1(interval_1 > 719) - 720; 
    integer_intervals_1 = [integer_intervals_1, interval_1];

    interval_6 = interval_1 + random_number_4; 
    interval_6(interval_6 > 719) = interval_6(interval_6 > 719) - 720;
    integer_intervals_6 = [integer_intervals_6, interval_6];
end

for i = 1:length(gg_0)
    interval_2 = start_points2(i) :end_points2(i) ;
    interval_2(interval_2 > 719) = interval_2(interval_2 > 719) - 720; 
    integer_intervals_2 = [integer_intervals_2, interval_2];

    interval_7 = interval_2 + random_number_5; 
    interval_7(interval_7 > 719) = interval_7(interval_7 > 719) - 720;
    integer_intervals_7 = [integer_intervals_7, interval_7];
end

for i = 1:length(ff_0)
    interval_1 = start_points1(i) :end_points1(i) ;
    interval_1(interval_1 > 719) = interval_1(interval_1 > 719) - 720; 
    integer_intervals_1 = [integer_intervals_1, interval_1];

    interval_8 = interval_1 + n; 
    interval_8(interval_8 > 719) = interval_8(interval_8 > 719) - 720;
    integer_intervals_8 = [integer_intervals_8, interval_8];
end


% 두 행렬을 합쳐서 하나의 행렬로 만듭니다.
sum_intervals = ...
    [integer_intervals_1, ...
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
    result_6 = [result_6; n, difference];
end

% result 변수에서 2열이 최소인 행을 찾습니다.
rows_with_min_6 = result_6(result_6(:, 2) == min(result_6(:, 2)), :);

k_6 = rows_with_min_6(:,1);


% 랜덤 인덱스 선택
random_index_6 = randi([1, length(k_6)]);

% 선택된 숫자
random_number_6 = k_6(random_index_6);

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
    interval_1 = start_points1(i) :end_points1(i) ;
    interval_1(interval_1 > 719) = interval_1(interval_1 > 719) - 720; 
    integer_intervals_1 = [integer_intervals_1, interval_1];
end


for i = 1:length(gg_0)
    interval_2 = start_points2(i) :end_points2(i) ;
    interval_2(interval_2 > 719) = interval_2(interval_2 > 719) - 720; 
    integer_intervals_2 = [integer_intervals_2, interval_2];

    interval_3 = interval_2 + random_number_1; 
    interval_3(interval_3 > 719) = interval_3(interval_3 > 719) - 720;
    integer_intervals_3 = [integer_intervals_3, interval_3];
end

for i = 1:length(ff_0)
    interval_1 = start_points1(i) :end_points1(i) ;
    interval_1(interval_1 > 719) = interval_1(interval_1 > 719) - 720; 
    integer_intervals_1 = [integer_intervals_1, interval_1];

    interval_4 = interval_1 + random_number_2; 
    interval_4(interval_4 > 719) = interval_4(interval_4 > 719) - 720;
    integer_intervals_4 = [integer_intervals_4, interval_4];
end

for i = 1:length(gg_0)
    interval_2 = start_points2(i) :end_points2(i) ;
    interval_2(interval_2 > 719) = interval_2(interval_2 > 719) - 720; 
    integer_intervals_2 = [integer_intervals_2, interval_2];

    interval_5 = interval_2 + random_number_3; 
    interval_5(interval_5 > 719) = interval_5(interval_5 > 719) - 720;
    integer_intervals_5 = [integer_intervals_5, interval_5];
end

for i = 1:length(ff_0)
    interval_1 = start_points1(i) :end_points1(i) ;
    interval_1(interval_1 > 719) = interval_1(interval_1 > 719) - 720; 
    integer_intervals_1 = [integer_intervals_1, interval_1];

    interval_6 = interval_1 + random_number_4; 
    interval_6(interval_6 > 719) = interval_6(interval_6 > 719) - 720;
    integer_intervals_6 = [integer_intervals_6, interval_6];
end

for i = 1:length(gg_0)
    interval_2 = start_points2(i) :end_points2(i) ;
    interval_2(interval_2 > 719) = interval_2(interval_2 > 719) - 720; 
    integer_intervals_2 = [integer_intervals_2, interval_2];

    interval_7 = interval_2 + random_number_5; 
    interval_7(interval_7 > 719) = interval_7(interval_7 > 719) - 720;
    integer_intervals_7 = [integer_intervals_7, interval_7];
end

for i = 1:length(ff_0)
    interval_1 = start_points1(i) :end_points1(i) ;
    interval_1(interval_1 > 719) = interval_1(interval_1 > 719) - 720; 
    integer_intervals_1 = [integer_intervals_1, interval_1];

    interval_8 = interval_1 + random_number_6; 
    interval_8(interval_8 > 719) = interval_8(interval_8 > 719) - 720;
    integer_intervals_8 = [integer_intervals_8, interval_8];
end


% 두 행렬을 합쳐서 하나의 행렬로 만듭니다.
sum_intervals = ...
    [integer_intervals_1, ...
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
    result_6 = [result_6; n, difference];
end

A_6 = length(unique_intervals) / 720 * 100;



%% 고궤도 추가 (4L / 4H)

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
    interval_1 = start_points1(i) :end_points1(i) ;
    interval_1(interval_1 > 719) = interval_1(interval_1 > 719) - 720; 
    integer_intervals_1 = [integer_intervals_1, interval_1];
end


for i = 1:length(gg_0)
    interval_2 = start_points2(i) :end_points2(i) ;
    interval_2(interval_2 > 719) = interval_2(interval_2 > 719) - 720; 
    integer_intervals_2 = [integer_intervals_2, interval_2];

    interval_3 = interval_2 + random_number_1; 
    interval_3(interval_3 > 719) = interval_3(interval_3 > 719) - 720;
    integer_intervals_3 = [integer_intervals_3, interval_3];
end

for i = 1:length(ff_0)
    interval_1 = start_points1(i) :end_points1(i) ;
    interval_1(interval_1 > 719) = interval_1(interval_1 > 719) - 720; 
    integer_intervals_1 = [integer_intervals_1, interval_1];

    interval_4 = interval_1 + random_number_2; 
    interval_4(interval_4 > 719) = interval_4(interval_4 > 719) - 720;
    integer_intervals_4 = [integer_intervals_4, interval_4];
end

for i = 1:length(gg_0)
    interval_2 = start_points2(i) :end_points2(i) ;
    interval_2(interval_2 > 719) = interval_2(interval_2 > 719) - 720; 
    integer_intervals_2 = [integer_intervals_2, interval_2];

    interval_5 = interval_2 + random_number_3; 
    interval_5(interval_5 > 719) = interval_5(interval_5 > 719) - 720;
    integer_intervals_5 = [integer_intervals_5, interval_5];
end

for i = 1:length(ff_0)
    interval_1 = start_points1(i) :end_points1(i) ;
    interval_1(interval_1 > 719) = interval_1(interval_1 > 719) - 720; 
    integer_intervals_1 = [integer_intervals_1, interval_1];

    interval_6 = interval_1 + random_number_4; 
    interval_6(interval_6 > 719) = interval_6(interval_6 > 719) - 720;
    integer_intervals_6 = [integer_intervals_6, interval_6];
end

for i = 1:length(gg_0)
    interval_2 = start_points2(i) :end_points2(i) ;
    interval_2(interval_2 > 719) = interval_2(interval_2 > 719) - 720; 
    integer_intervals_2 = [integer_intervals_2, interval_2];

    interval_7 = interval_2 + random_number_5; 
    interval_7(interval_7 > 719) = interval_7(interval_7 > 719) - 720;
    integer_intervals_7 = [integer_intervals_7, interval_7];
end

for i = 1:length(ff_0)
    interval_1 = start_points1(i) :end_points1(i) ;
    interval_1(interval_1 > 719) = interval_1(interval_1 > 719) - 720; 
    integer_intervals_1 = [integer_intervals_1, interval_1];

    interval_8 = interval_1 + random_number_6; 
    interval_8(interval_8 > 719) = interval_8(interval_8 > 719) - 720;
    integer_intervals_8 = [integer_intervals_8, interval_8];
end

for i = 1:length(gg_0)
    interval_2 = start_points2(i) :end_points2(i) ;
    interval_2(interval_2 > 719) = interval_2(interval_2 > 719) - 720; 
    integer_intervals_2 = [integer_intervals_2, interval_2];

    interval_9 = interval_2 + n; 
    interval_9(interval_9 > 719) = interval_9(interval_9 > 719) - 720;
    integer_intervals_9 = [integer_intervals_9, interval_9];
end


% 두 행렬을 합쳐서 하나의 행렬로 만듭니다.
sum_intervals = ...
    [integer_intervals_1, ...
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
    result_7 = [result_7; n, difference];
end

% result 변수에서 2열이 최소인 행을 찾습니다.
rows_with_min_7 = result_7(result_7(:, 2) == min(result_7(:, 2)), :);

k_7 = rows_with_min_7(:,1);


% 랜덤 인덱스 선택
random_index_7 = randi([1, length(k_7)]);

% 선택된 숫자
random_number_7 = k_7(random_index_7);

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
    interval_1 = start_points1(i) :end_points1(i) ;
    interval_1(interval_1 > 719) = interval_1(interval_1 > 719) - 720; 
    integer_intervals_1 = [integer_intervals_1, interval_1];
end


for i = 1:length(gg_0)
    interval_2 = start_points2(i) :end_points2(i) ;
    interval_2(interval_2 > 719) = interval_2(interval_2 > 719) - 720; 
    integer_intervals_2 = [integer_intervals_2, interval_2];

    interval_3 = interval_2 + random_number_1; 
    interval_3(interval_3 > 719) = interval_3(interval_3 > 719) - 720;
    integer_intervals_3 = [integer_intervals_3, interval_3];
end

for i = 1:length(ff_0)
    interval_1 = start_points1(i) :end_points1(i) ;
    interval_1(interval_1 > 719) = interval_1(interval_1 > 719) - 720; 
    integer_intervals_1 = [integer_intervals_1, interval_1];

    interval_4 = interval_1 + random_number_2; 
    interval_4(interval_4 > 719) = interval_4(interval_4 > 719) - 720;
    integer_intervals_4 = [integer_intervals_4, interval_4];
end

for i = 1:length(gg_0)
    interval_2 = start_points2(i) :end_points2(i) ;
    interval_2(interval_2 > 719) = interval_2(interval_2 > 719) - 720; 
    integer_intervals_2 = [integer_intervals_2, interval_2];

    interval_5 = interval_2 + random_number_3; 
    interval_5(interval_5 > 719) = interval_5(interval_5 > 719) - 720;
    integer_intervals_5 = [integer_intervals_5, interval_5];
end

for i = 1:length(ff_0)
    interval_1 = start_points1(i) :end_points1(i) ;
    interval_1(interval_1 > 719) = interval_1(interval_1 > 719) - 720; 
    integer_intervals_1 = [integer_intervals_1, interval_1];

    interval_6 = interval_1 + random_number_4; 
    interval_6(interval_6 > 719) = interval_6(interval_6 > 719) - 720;
    integer_intervals_6 = [integer_intervals_6, interval_6];
end

for i = 1:length(gg_0)
    interval_2 = start_points2(i) :end_points2(i) ;
    interval_2(interval_2 > 719) = interval_2(interval_2 > 719) - 720; 
    integer_intervals_2 = [integer_intervals_2, interval_2];

    interval_7 = interval_2 + random_number_5; 
    interval_7(interval_7 > 719) = interval_7(interval_7 > 719) - 720;
    integer_intervals_7 = [integer_intervals_7, interval_7];
end

for i = 1:length(ff_0)
    interval_1 = start_points1(i) :end_points1(i) ;
    interval_1(interval_1 > 719) = interval_1(interval_1 > 719) - 720; 
    integer_intervals_1 = [integer_intervals_1, interval_1];

    interval_8 = interval_1 + random_number_6; 
    interval_8(interval_8 > 719) = interval_8(interval_8 > 719) - 720;
    integer_intervals_8 = [integer_intervals_8, interval_8];
end

for i = 1:length(gg_0)
    interval_2 = start_points2(i) :end_points2(i) ;
    interval_2(interval_2 > 719) = interval_2(interval_2 > 719) - 720; 
    integer_intervals_2 = [integer_intervals_2, interval_2];

    interval_9 = interval_2 + random_number_7; 
    interval_9(interval_9 > 719) = interval_9(interval_9 > 719) - 720;
    integer_intervals_9 = [integer_intervals_9, interval_9];
end


% 두 행렬을 합쳐서 하나의 행렬로 만듭니다.
sum_intervals = ...
    [integer_intervals_1, ...
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
    result_7 = [result_7; n, difference];
end

A_7 = length(unique_intervals) / 720 * 100;


%% 고궤도 추가 (4L / 5H)

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
    interval_1 = start_points1(i) :end_points1(i) ;
    interval_1(interval_1 > 719) = interval_1(interval_1 > 719) - 720; 
    integer_intervals_1 = [integer_intervals_1, interval_1];
end


for i = 1:length(gg_0)
    interval_2 = start_points2(i) :end_points2(i) ;
    interval_2(interval_2 > 719) = interval_2(interval_2 > 719) - 720; 
    integer_intervals_2 = [integer_intervals_2, interval_2];

    interval_3 = interval_2 + random_number_1; 
    interval_3(interval_3 > 719) = interval_3(interval_3 > 719) - 720;
    integer_intervals_3 = [integer_intervals_3, interval_3];
end

for i = 1:length(ff_0)
    interval_1 = start_points1(i) :end_points1(i) ;
    interval_1(interval_1 > 719) = interval_1(interval_1 > 719) - 720; 
    integer_intervals_1 = [integer_intervals_1, interval_1];

    interval_4 = interval_1 + random_number_2; 
    interval_4(interval_4 > 719) = interval_4(interval_4 > 719) - 720;
    integer_intervals_4 = [integer_intervals_4, interval_4];
end

for i = 1:length(gg_0)
    interval_2 = start_points2(i) :end_points2(i) ;
    interval_2(interval_2 > 719) = interval_2(interval_2 > 719) - 720; 
    integer_intervals_2 = [integer_intervals_2, interval_2];

    interval_5 = interval_2 + random_number_3; 
    interval_5(interval_5 > 719) = interval_5(interval_5 > 719) - 720;
    integer_intervals_5 = [integer_intervals_5, interval_5];
end

for i = 1:length(ff_0)
    interval_1 = start_points1(i) :end_points1(i) ;
    interval_1(interval_1 > 719) = interval_1(interval_1 > 719) - 720; 
    integer_intervals_1 = [integer_intervals_1, interval_1];

    interval_6 = interval_1 + random_number_4; 
    interval_6(interval_6 > 719) = interval_6(interval_6 > 719) - 720;
    integer_intervals_6 = [integer_intervals_6, interval_6];
end

for i = 1:length(gg_0)
    interval_2 = start_points2(i) :end_points2(i) ;
    interval_2(interval_2 > 719) = interval_2(interval_2 > 719) - 720; 
    integer_intervals_2 = [integer_intervals_2, interval_2];

    interval_7 = interval_2 + random_number_5; 
    interval_7(interval_7 > 719) = interval_7(interval_7 > 719) - 720;
    integer_intervals_7 = [integer_intervals_7, interval_7];
end

for i = 1:length(ff_0)
    interval_1 = start_points1(i) :end_points1(i) ;
    interval_1(interval_1 > 719) = interval_1(interval_1 > 719) - 720; 
    integer_intervals_1 = [integer_intervals_1, interval_1];

    interval_8 = interval_1 + random_number_6; 
    interval_8(interval_8 > 719) = interval_8(interval_8 > 719) - 720;
    integer_intervals_8 = [integer_intervals_8, interval_8];
end

for i = 1:length(gg_0)
    interval_2 = start_points2(i) :end_points2(i) ;
    interval_2(interval_2 > 719) = interval_2(interval_2 > 719) - 720; 
    integer_intervals_2 = [integer_intervals_2, interval_2];

    interval_9 = interval_2 + random_number_7; 
    interval_9(interval_9 > 719) = interval_9(interval_9 > 719) - 720;
    integer_intervals_9 = [integer_intervals_9, interval_9];
end

for i = 1:length(gg_0)
    interval_2 = start_points2(i) :end_points2(i) ;
    interval_2(interval_2 > 719) = interval_2(interval_2 > 719) - 720; 
    integer_intervals_2 = [integer_intervals_2, interval_2];

    interval_10 = interval_2 + n; 
    interval_10(interval_10 > 719) = interval_10(interval_10 > 719) - 720;
    integer_intervals_10 = [integer_intervals_10, interval_10];
end


% 두 행렬을 합쳐서 하나의 행렬로 만듭니다.
sum_intervals = ...
    [integer_intervals_1, ...
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
    result_8 = [result_8; n, difference];
end

% result 변수에서 2열이 최소인 행을 찾습니다.
rows_with_min_8 = result_8(result_8(:, 2) == min(result_8(:, 2)), :);

k_8 = rows_with_min_8(:,1);


% 랜덤 인덱스 선택
random_index_8 = randi([1, length(k_8)]);

% 선택된 숫자
random_number_8 = k_8(random_index_8);

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
    interval_1 = start_points1(i) :end_points1(i) ;
    interval_1(interval_1 > 719) = interval_1(interval_1 > 719) - 720; 
    integer_intervals_1 = [integer_intervals_1, interval_1];
end


for i = 1:length(gg_0)
    interval_2 = start_points2(i) :end_points2(i) ;
    interval_2(interval_2 > 719) = interval_2(interval_2 > 719) - 720; 
    integer_intervals_2 = [integer_intervals_2, interval_2];

    interval_3 = interval_2 + random_number_1; 
    interval_3(interval_3 > 719) = interval_3(interval_3 > 719) - 720;
    integer_intervals_3 = [integer_intervals_3, interval_3];
end


for i = 1:length(ff_0)
    interval_1 = start_points1(i) :end_points1(i) ;
    interval_1(interval_1 > 719) = interval_1(interval_1 > 719) - 720; 
    integer_intervals_1 = [integer_intervals_1, interval_1];

    interval_4 = interval_1 + random_number_2; 
    interval_4(interval_4 > 719) = interval_4(interval_4 > 719) - 720;
    integer_intervals_4 = [integer_intervals_4, interval_4];
end


for i = 1:length(gg_0)
    interval_2 = start_points2(i) :end_points2(i) ;
    interval_2(interval_2 > 719) = interval_2(interval_2 > 719) - 720; 
    integer_intervals_2 = [integer_intervals_2, interval_2];

    interval_5 = interval_2 + random_number_3; 
    interval_5(interval_5 > 719) = interval_5(interval_5 > 719) - 720;
    integer_intervals_5 = [integer_intervals_5, interval_5];
end


for i = 1:length(ff_0)
    interval_1 = start_points1(i) :end_points1(i) ;
    interval_1(interval_1 > 719) = interval_1(interval_1 > 719) - 720; 
    integer_intervals_1 = [integer_intervals_1, interval_1];

    interval_6 = interval_1 + random_number_4; 
    interval_6(interval_6 > 719) = interval_6(interval_6 > 719) - 720;
    integer_intervals_6 = [integer_intervals_6, interval_6];
end


for i = 1:length(gg_0)
    interval_2 = start_points2(i) :end_points2(i) ;
    interval_2(interval_2 > 719) = interval_2(interval_2 > 719) - 720; 
    integer_intervals_2 = [integer_intervals_2, interval_2];

    interval_7 = interval_2 + random_number_5; 
    interval_7(interval_7 > 719) = interval_7(interval_7 > 719) - 720;
    integer_intervals_7 = [integer_intervals_7, interval_7];
end


for i = 1:length(ff_0)
    interval_1 = start_points1(i) :end_points1(i) ;
    interval_1(interval_1 > 719) = interval_1(interval_1 > 719) - 720; 
    integer_intervals_1 = [integer_intervals_1, interval_1];

    interval_8 = interval_1 + random_number_6; 
    interval_8(interval_8 > 719) = interval_8(interval_8 > 719) - 720;
    integer_intervals_8 = [integer_intervals_8, interval_8];
end


for i = 1:length(gg_0)
    interval_2 = start_points2(i) :end_points2(i) ;
    interval_2(interval_2 > 719) = interval_2(interval_2 > 719) - 720; 
    integer_intervals_2 = [integer_intervals_2, interval_2];

    interval_9 = interval_2 + random_number_7; 
    interval_9(interval_9 > 719) = interval_9(interval_9 > 719) - 720;
    integer_intervals_9 = [integer_intervals_9, interval_9];
end


for i = 1:length(gg_0)
    interval_2 = start_points2(i) :end_points2(i) ;
    interval_2(interval_2 > 719) = interval_2(interval_2 > 719) - 720; 
    integer_intervals_2 = [integer_intervals_2, interval_2];

    interval_10 = interval_2 + random_number_8; 
    interval_10(interval_10 > 719) = interval_10(interval_10 > 719) - 720;
    integer_intervals_10 = [integer_intervals_10, interval_10];
end


% 두 행렬을 합쳐서 하나의 행렬로 만듭니다.
sum_intervals = ...
    [integer_intervals_1, ...
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
    result_8 = [result_8; n, difference];
end

A_8 = length(unique_intervals) / 720 * 100


%% 저궤도 추가 (5L / 5H)

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
    interval_1 = start_points1(i) :end_points1(i) ;
    interval_1(interval_1 > 719) = interval_1(interval_1 > 719) - 720; 
    integer_intervals_1 = [integer_intervals_1, interval_1];
end


for i = 1:length(gg_0)
    interval_2 = start_points2(i) :end_points2(i) ;
    interval_2(interval_2 > 719) = interval_2(interval_2 > 719) - 720; 
    integer_intervals_2 = [integer_intervals_2, interval_2];

    interval_3 = interval_2 + random_number_1; 
    interval_3(interval_3 > 719) = interval_3(interval_3 > 719) - 720;
    integer_intervals_3 = [integer_intervals_3, interval_3];
end

for i = 1:length(ff_0)
    interval_1 = start_points1(i) :end_points1(i) ;
    interval_1(interval_1 > 719) = interval_1(interval_1 > 719) - 720; 
    integer_intervals_1 = [integer_intervals_1, interval_1];

    interval_4 = interval_1 + random_number_2; 
    interval_4(interval_4 > 719) = interval_4(interval_4 > 719) - 720;
    integer_intervals_4 = [integer_intervals_4, interval_4];
end

for i = 1:length(gg_0)
    interval_2 = start_points2(i) :end_points2(i) ;
    interval_2(interval_2 > 719) = interval_2(interval_2 > 719) - 720; 
    integer_intervals_2 = [integer_intervals_2, interval_2];

    interval_5 = interval_2 + random_number_3; 
    interval_5(interval_5 > 719) = interval_5(interval_5 > 719) - 720;
    integer_intervals_5 = [integer_intervals_5, interval_5];
end

for i = 1:length(ff_0)
    interval_1 = start_points1(i) :end_points1(i) ;
    interval_1(interval_1 > 719) = interval_1(interval_1 > 719) - 720; 
    integer_intervals_1 = [integer_intervals_1, interval_1];

    interval_6 = interval_1 + random_number_4; 
    interval_6(interval_6 > 719) = interval_6(interval_6 > 719) - 720;
    integer_intervals_6 = [integer_intervals_6, interval_6];
end

for i = 1:length(gg_0)
    interval_2 = start_points2(i) :end_points2(i) ;
    interval_2(interval_2 > 719) = interval_2(interval_2 > 719) - 720; 
    integer_intervals_2 = [integer_intervals_2, interval_2];

    interval_7 = interval_2 + random_number_5; 
    interval_7(interval_7 > 719) = interval_7(interval_7 > 719) - 720;
    integer_intervals_7 = [integer_intervals_7, interval_7];
end

for i = 1:length(ff_0)
    interval_1 = start_points1(i) :end_points1(i) ;
    interval_1(interval_1 > 719) = interval_1(interval_1 > 719) - 720; 
    integer_intervals_1 = [integer_intervals_1, interval_1];

    interval_8 = interval_1 + random_number_6; 
    interval_8(interval_8 > 719) = interval_8(interval_8 > 719) - 720;
    integer_intervals_8 = [integer_intervals_8, interval_8];
end

for i = 1:length(gg_0)
    interval_2 = start_points2(i) :end_points2(i) ;
    interval_2(interval_2 > 719) = interval_2(interval_2 > 719) - 720; 
    integer_intervals_2 = [integer_intervals_2, interval_2];

    interval_9 = interval_2 + random_number_7; 
    interval_9(interval_9 > 719) = interval_9(interval_9 > 719) - 720;
    integer_intervals_9 = [integer_intervals_9, interval_9];
end

for i = 1:length(gg_0)
    interval_2 = start_points2(i) :end_points2(i) ;
    interval_2(interval_2 > 719) = interval_2(interval_2 > 719) - 720; 
    integer_intervals_2 = [integer_intervals_2, interval_2];

    interval_10 = interval_2 + random_number_8; 
    interval_10(interval_10 > 719) = interval_10(interval_10 > 719) - 720;
    integer_intervals_10 = [integer_intervals_10, interval_10];
end

for i = 1:length(gg_0)
    interval_2 = start_points2(i) :end_points2(i) ;
    interval_2(interval_2 > 719) = interval_2(interval_2 > 719) - 720; 
    integer_intervals_2 = [integer_intervals_2, interval_2];

    interval_11 = interval_2 + n; 
    interval_11(interval_11 > 719) = interval_11(interval_11 > 719) - 720;
    integer_intervals_11 = [integer_intervals_11, interval_11];
end


% 두 행렬을 합쳐서 하나의 행렬로 만듭니다.
sum_intervals = ...
    [integer_intervals_1, ...
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
    result_9 = [result_9; n, difference];
end

% result 변수에서 2열이 최소인 행을 찾습니다.
rows_with_min_9 = result_9(result_9(:, 2) == min(result_9(:, 2)), :);

k_9 = rows_with_min_9(:,1);


% 랜덤 인덱스 선택
random_index_9 = randi([1, length(k_9)]);

% 선택된 숫자
random_number_9 = k_9(random_index_9);

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
    interval_1 = start_points1(i) :end_points1(i) ;
    interval_1(interval_1 > 719) = interval_1(interval_1 > 719) - 720; 
    integer_intervals_1 = [integer_intervals_1, interval_1];
end


for i = 1:length(gg_0)
    interval_2 = start_points2(i) :end_points2(i) ;
    interval_2(interval_2 > 719) = interval_2(interval_2 > 719) - 720; 
    integer_intervals_2 = [integer_intervals_2, interval_2];

    interval_3 = interval_2 + random_number_1; 
    interval_3(interval_3 > 719) = interval_3(interval_3 > 719) - 720;
    integer_intervals_3 = [integer_intervals_3, interval_3];
end


for i = 1:length(ff_0)
    interval_1 = start_points1(i) :end_points1(i) ;
    interval_1(interval_1 > 719) = interval_1(interval_1 > 719) - 720; 
    integer_intervals_1 = [integer_intervals_1, interval_1];

    interval_4 = interval_1 + random_number_2; 
    interval_4(interval_4 > 719) = interval_4(interval_4 > 719) - 720;
    integer_intervals_4 = [integer_intervals_4, interval_4];
end


for i = 1:length(gg_0)
    interval_2 = start_points2(i) :end_points2(i) ;
    interval_2(interval_2 > 719) = interval_2(interval_2 > 719) - 720; 
    integer_intervals_2 = [integer_intervals_2, interval_2];

    interval_5 = interval_2 + random_number_3; 
    interval_5(interval_5 > 719) = interval_5(interval_5 > 719) - 720;
    integer_intervals_5 = [integer_intervals_5, interval_5];
end


for i = 1:length(ff_0)
    interval_1 = start_points1(i) :end_points1(i) ;
    interval_1(interval_1 > 719) = interval_1(interval_1 > 719) - 720; 
    integer_intervals_1 = [integer_intervals_1, interval_1];

    interval_6 = interval_1 + random_number_4; 
    interval_6(interval_6 > 719) = interval_6(interval_6 > 719) - 720;
    integer_intervals_6 = [integer_intervals_6, interval_6];
end


for i = 1:length(gg_0)
    interval_2 = start_points2(i) :end_points2(i) ;
    interval_2(interval_2 > 719) = interval_2(interval_2 > 719) - 720; 
    integer_intervals_2 = [integer_intervals_2, interval_2];

    interval_7 = interval_2 + random_number_5; 
    interval_7(interval_7 > 719) = interval_7(interval_7 > 719) - 720;
    integer_intervals_7 = [integer_intervals_7, interval_7];
end


for i = 1:length(ff_0)
    interval_1 = start_points1(i) :end_points1(i) ;
    interval_1(interval_1 > 719) = interval_1(interval_1 > 719) - 720; 
    integer_intervals_1 = [integer_intervals_1, interval_1];

    interval_8 = interval_1 + random_number_6; 
    interval_8(interval_8 > 719) = interval_8(interval_8 > 719) - 720;
    integer_intervals_8 = [integer_intervals_8, interval_8];
end


for i = 1:length(gg_0)
    interval_2 = start_points2(i) :end_points2(i) ;
    interval_2(interval_2 > 719) = interval_2(interval_2 > 719) - 720; 
    integer_intervals_2 = [integer_intervals_2, interval_2];

    interval_9 = interval_2 + random_number_7; 
    interval_9(interval_9 > 719) = interval_9(interval_9 > 719) - 720;
    integer_intervals_9 = [integer_intervals_9, interval_9];
end


for i = 1:length(gg_0)
    interval_2 = start_points2(i) :end_points2(i) ;
    interval_2(interval_2 > 719) = interval_2(interval_2 > 719) - 720; 
    integer_intervals_2 = [integer_intervals_2, interval_2];

    interval_10 = interval_2 + random_number_8; 
    interval_10(interval_10 > 719) = interval_10(interval_10 > 719) - 720;
    integer_intervals_10 = [integer_intervals_10, interval_10];
end

for i = 1:length(gg_0)
    interval_2 = start_points2(i) :end_points2(i) ;
    interval_2(interval_2 > 719) = interval_2(interval_2 > 719) - 720; 
    integer_intervals_2 = [integer_intervals_2, interval_2];

    interval_11 = interval_2 + random_number_9; 
    interval_11(interval_11 > 719) = interval_11(interval_11 > 719) - 720;
    integer_intervals_11 = [integer_intervals_11, interval_11];
end

% 두 행렬을 합쳐서 하나의 행렬로 만듭니다.
sum_intervals = ...
    [integer_intervals_1, ...
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
    result_9 = [result_9; n, difference];
end

A_10 = length(unique_intervals) / 720 * 100

%%

random_N(1) = (random_number_1);
random_N(2) = (random_number_2);
random_N(3) = (random_number_3);
random_N(4) = (random_number_4);
random_N(5) = (random_number_5);
random_N(6) = (random_number_6);
random_N(7) = (random_number_7);
random_N(8) = (random_number_8)
% random_N(9) = (random_number_9)

  if A_8 >= 30


        break;
    end
 trials = trials + 1; 
end

% 시뮬레이션 종료 시간 기록
elapsed_time = toc(start_time);

% 결과 출력
fprintf('시행 횟수: %d\n', trials);
fprintf('걸린 시간(초): %.2f\n', elapsed_time);







%% 제궤도 그래프

n_1 = [0, random_number_2, random_number_4, random_number_6]+1;

%ff_0 = seedsatellite
ff_1 = ff_0 + n_1(1);
ff_2 = ff_0 + n_1(2);
ff_3 = ff_0 + n_1(3);
ff_4 = ff_0 + n_1(4);
% ff_5 = ff_0 + n_1(5);

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

coverage1 = length(unique_x_coords1) / L * 100;

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
RV = [RV; current_value, end_index - start_index+1];


%% 빠진 숫자 찾기

% 주어진 행렬
matrix = RV ;

% 현재 1열에서 빠진 숫자를 찾기
missing_numbers = setdiff(0:L-1, matrix(:,1));

% 새롭게 추가된 행에 대해 0으로 입력
new_rows = zeros(length(missing_numbers), 2);
new_rows(:, 1) = missing_numbers;

% 새로운 행렬 생성
new_matrix = [matrix; new_rows];


% 1열을 기준으로 행렬을 정렬
sorted_matrix1 = sortrows(new_matrix, 1);

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




%% 고궤도 그래프

n_2 = [random_number_1, random_number_3, random_number_5, random_number_7, random_number_8, random_number_9]+1;

%ff_0 = seedsatellite
gg_1 = gg_0 + n_2(1);
gg_2 = gg_0 + n_2(2);
gg_3 = gg_0 + n_2(3);
gg_4 = gg_0 + n_2(4);
gg_5 = gg_0 + n_2(5);
gg_6 = gg_0 + n_2(6);

% x축 범위 설정
x_range = 0:L-1;

for idx = 1:length(n_2)
    gg = eval(['gg_', num2str(idx)]);
    x_coords2{idx} = []; % 해당 작업공간의 x좌표를 저장할 배열 초기화

    for i = 1:size(gg,1)
        start_time = gg(i, 1);
        end_time = gg(i, 2);

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
          x_coords2{idx} = [x_coords2{idx}, x_range(start_index:L),(0:end_index)];

        else
          x_coords2{idx} = [x_coords2{idx}, x_range(start_index:end_index)];
        end



    end               
end

% 모든 작업공간의 x 좌표를 모은 배열 초기화
all_x_coords2 = [];

%% 각 작업공간의 x 좌표를 all_x_coords1에 추가하고 작은 수부터 오름차순으로 정렬
for idx = 1:length(n_2)
    all_x_coords2 = [all_x_coords2, x_coords2{idx}];
end

all_x_coords2 = sort(all_x_coords2);

% 중복된 값을 제거하고 유일한 값들만 추출
unique_x_coords2 = unique(all_x_coords2);

coverage2 = length(unique_x_coords2) / L * 100;

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
RV = [RV; current_value, end_index - start_index+1];


%% 빠진 숫자 찾기

% 주어진 행렬
matrix = RV ;

% 현재 1열에서 빠진 숫자를 찾기
missing_numbers = setdiff(0:L-1, matrix(:,1));

% 새롭게 추가된 행에 대해 0으로 입력
new_rows = zeros(length(missing_numbers), 2);
new_rows(:, 1) = missing_numbers;

% 새로운 행렬 생성
new_matrix = [matrix; new_rows];


% 1열을 기준으로 행렬을 정렬
sorted_matrix2 = sortrows(new_matrix, 1);

% Acess Profile

% x축 범위 설정
x_range = 0:L-1;

figure;
subplot(3, 1, 1);

hex_color2 = '#cc4e01';

red2 = hex2dec(hex_color2(2:3)) / 255; % R 값
green2 = hex2dec(hex_color2(4:5)) / 255; % G 값
blue2 = hex2dec(hex_color2(6:7)) / 255; % B 값

hex_color22 = '#fe9882';

red22 = hex2dec(hex_color22(2:3)) / 255; % R 값
green22 = hex2dec(hex_color22(4:5)) / 255; % G 값
blue22 = hex2dec(hex_color22(6:7)) / 255; % B 값


% 각 구간별로 y=1인 위치 표시
for i = 1:size(gg_0, 1)
    if neg2(i) < 0
        start_time_new = L + neg2(i,1);
        end_time_new = L + neg2(i,2);
    else 
        start_time_new = gg_0(i, 1);
        end_time_new = gg_0(i, 2);
    end
    start_time = gg_0(i, 1);
    end_time = gg_0(i, 2);

    if end_time_new > L-1
        end_time_new = L-1;
    end

    % 시작 시간과 종료 시간에 해당하는 x 인덱스 찾기
    start_index = find(x_range >= start_time, 1);
    end_index = find(x_range <= end_time, 1, 'last');
    stn = start_time_new;
    etn = end_time_new;

    % y 값 설정
    y_values2a = zeros(size(x_range));
    y_values2a(start_index:end_index) = 1;
    if neg2(i) < 0
        y_values2a(stn:etn) = 1;
    end

    % 그래프 그리기
    plot(x_range, y_values2a, '-o r', 'LineWidth', 1, 'MarkerSize', 1, 'Color', [red2, green2, blue2]);
    hold on;
    area(x_range, y_values2a, 'FaceColor', [red22, green22, blue22]);
    plot(x_range, y_values2a, 'Color', [red2, green2, blue2], 'LineWidth', 1.5); 
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
    y_values2p = zeros(size(x_range));
    y_values2p(n_2(1,:)) = 1;


plot(x_range, y_values2p, 'Color', [red2, green2, blue2], 'LineWidth', 1, 'MarkerSize', 1);
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
plot(sorted_matrix2(:,1), sorted_matrix2(:,2), 'Color', [red2, green2, blue2], 'LineWidth', 1, 'MarkerSize', 1);
title('Coverage Timeline');
hold on;
area(sorted_matrix2(:,1), sorted_matrix2(:,2), 'FaceColor', [red22, green22, blue22]);
plot(sorted_matrix2(:,1), sorted_matrix2(:,2), 'Color', [red2, green2, blue2], 'LineWidth', 1);

xlim([0, L-1]);
ylim([0, 3]);

% 그래프 세부 설정
xlabel('n');
ylabel('bj[n]');
grid on;
xlim([0, L-1]);


%%

% 새로운 그래프에 두 개의 플롯을 함께 보기
figure3 = figure;
subplot(3, 1, 1);

hold on;


% ff
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
    hold on;
    plot(x_range, y_values1a, '-o r', 'LineWidth', 1, 'MarkerSize', 1, 'Color', [red1, green1, blue1]);
    area(x_range, y_values1a, 'FaceColor', [red11, green11, blue11], 'FaceAlpha', 0.7);
    plot(x_range, y_values1a, 'Color', [red1, green1, blue1], 'LineWidth', 1.5); 
end

% gg
for i = 1:size(gg_0, 1)
    if neg2(i) < 0
        start_time_new = L + neg2(i,1);
        end_time_new = L + neg2(i,2);
    else 
        start_time_new = gg_0(i, 1);
        end_time_new = gg_0(i, 1);
    end
    start_time = gg_0(i, 1);
    end_time = gg_0(i, 2);

    if end_time_new > L-1
        end_time_new = L-1;
    end

    % 시작 시간과 종료 시간에 해당하는 x 인덱스 찾기
    start_index = find(x_range >= start_time, 1);
    end_index = find(x_range <= end_time, 1, 'last');
    stn = start_time_new;
    etn = end_time_new;

    % y 값 설정
    y_values2a = zeros(size(x_range));
    y_values2a(start_index:end_index) = 1;
    if neg2(i) < 0
        y_values2a(stn:etn) = 1;
    end

    % 그래프 그리기
    plot(x_range, y_values2a, '-o r', 'LineWidth', 1, 'MarkerSize', 1, 'Color', [red2, green2, blue2]);
    hold on;
    area(x_range, y_values2a, 'FaceColor', [red22, green22, blue22], 'FaceAlpha', 0.7);
    plot(x_range, y_values2a, 'Color', [red2, green2, blue2], 'LineWidth', 1.5);
end


hold on; % 여러 플롯을 함께 그리기 위해 hold on 사용
plot(x_range, y_values1a, 'Color', [red1, green1, blue1]); % 첫 번째 플롯 복사
plot(x_range, y_values2a, 'Color', [red2, green2, blue2]); % 두 번째 플롯 복사
xlim([0, L-1]);
ylim([0, 1]);
yticks([0, 1]);

% 그래프 세부 설정
ylabel('vo,j[n]');
title('Seed Satellite Access Profile');
grid on;



subplot(3, 1, 2);
hold on; % 여러 플롯을 함께 그리기 위해 hold on 사용
plot(x_range, y_values1p, 'Color', [red1, green1, blue1]); % 첫 번째 플롯 복사
plot(x_range, y_values2p, 'Color', [red2, green2, blue2]); % 두 번째 플롯 복사
xlim([0, L-1]);
ylim([0, 1]);
yticks([0, 1]);

% 그래프 세부 설정
ylabel('x[n]');
title('Constellation Pattern Vector');
grid on;



subplot(3, 1, 3);
hold on; % 여러 플롯을 함께 그리기 위해 hold on 사용
plot(sorted_matrix1(:,1), sorted_matrix1(:,2), 'Color', [red1, green1, blue1]); % 첫 번째 플롯 복사
hold on;
area(sorted_matrix1(:,1), sorted_matrix1(:,2), 'FaceColor', [red11, green11, blue11], 'FaceAlpha', 0.7);
plot(sorted_matrix1(:,1), sorted_matrix1(:,2), 'Color', [red1, green1, blue1], 'LineWidth', 1.5);

plot(sorted_matrix2(:,1), sorted_matrix2(:,2), 'Color', [red2, green2, blue2]); % 두 번째 플롯 복사
hold on;
area(sorted_matrix2(:,1), sorted_matrix2(:,2), 'FaceColor', [red22, green22, blue22], 'FaceAlpha', 0.7);
plot(sorted_matrix2(:,1), sorted_matrix2(:,2), 'Color', [red2, green2, blue2], 'LineWidth', 1.5);

xlim([0, L-1]);
ylim([0, 3]);

% 그래프 세부 설정
ylabel('n');
ylabel('bj[n]');
title('Coverage Timeline');
grid on;

hold off; % hold 해제


% 새로운 행렬 생성 (첫 번째 열은 아무 변경 없음)
sorted_matrix3 = zeros(720, 2);  % 결과 행렬 초기화
sorted_matrix3(:, 1) = sorted_matrix1(:, 1);  % 1열은 첫 번째 행렬에서 복사

% 2열의 OR 연산을 사용하여 값 결정
sorted_matrix3(:, 2) = sorted_matrix1(:, 2) | sorted_matrix2(:, 2);

sorted_matrix3(:,2) = sorted_matrix1(:,2) + sorted_matrix2(:,2);

count1 = sum(sorted_matrix1(:, 2) >= 1); 
count2 = sum(sorted_matrix2(:, 2) >= 1); 
count = sum(sorted_matrix3(:, 2) >= 1); 

coverage_total = (count) / L * 100

all_x_coords3 = unique(sort(all_x_coords1(:,2) + all_x_coords2(:,2)));

