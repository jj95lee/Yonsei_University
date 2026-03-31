clc
clear all

L= 717;
shift = 0;

% ff_0 -> n 으로 변환
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

% ff_0 -> n 으로 변환
ff_0 = floor(ff_00/(86024/L)) - shift;

% 구간의 시작점과 끝점을 변수에 저장
start_points = ff_0(:, 1);
end_points = ff_0(:, 2);


% 어떤 수를 골랐을 때 중복이 몇인지
result = [];
result_1 = [];
result_2 = [];
result_3 = [];

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
sum_intervals = [integer_intervals_1, integer_intervals_2];

% 행렬을 오름차순으로 정렬하고 중복된 요소를 제거합니다.
unique_intervals = unique(sort(sum_intervals));

% 두 행렬의 길이 차이 계산
difference = length(sum_intervals) - length(unique_intervals) ;

% 결과를 작업 공간에 저장
result_1 = [result_1; n, difference];

end

A_2 = length(unique_intervals) / 717 * 100;

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


% % 작업공간에 저장
% RAAN(k+1) = rightAscensionOfAscendingNode;   
% Ta(k+1) = trueAnomaly;
% n_pos(k+1) = round(720/ N *k,0); 

% 작업공간에 저장
% K_1(j) = length(n);   


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

A_3 = length(unique_intervals) / 717 * 100;

% result 변수에서 2열이 최소인 행을 찾습니다.
rows_with_min_2 = result_2(result_2(:, 2) == min(result_2(:, 2)), :);

k_2 = rows_with_min_2(:,1);


% result = [result; k_1(j), n];








