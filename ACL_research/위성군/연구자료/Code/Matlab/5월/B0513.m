clear all

trials = 0;
start_time = tic;
L = 720;

while true

% korea 
ff_00 = ...
[...
32032.79297	34606.58203
43446.85547	46222.20703
55434.55078	58151.13281
67135.60547	69899.35547
78948.39844	80574.08203
...		
];

gg_00 = ...
[...
8548.886719	13105.89844
24845.97656	29186.30859
78420.52734	82803.10547
...		
];

% example
hh_00 = ...
[...
26797.85156	28499.12109
37473.22266	40238.67188
49309.21875	51894.72656
61331.30859	64021.81641
72886.40625	75528.33984
...		
];

ii_00 = ...
[...
65.9765625	4605
16671.26953	21215.91797
32836.46484	36689.35547
70896.62109	74453.84766
86094.66797	86400

...		
];


% 저궤도

% spot1 korea
% 31963.65234	34335.9375
% 43415.21484	46260.05859
% 55296.21094	58167.65625
% 67180.66406	69859.10156
% 79699.33594	80474.94141

% spot 2 example
% 31883.02734	34448.20313
% 43490.97656	46362.36328
% 55453.76953	58329.78516
% 67360.3125	70037.69531


% 고궤도

% spot 1 korea
% 27825.41016	31367.40234
% 43297.20703	47883.69141
% 59740.83984	64309.98047
% 76308.45703	79650.87891

% spot 2 example
% 27642.42188	31526.66016
% 43402.08984	48035.74219
% 59971.40625	64554.25781
% 76648.125	79865.56641


% ff_0 -> n 으로 변환
ff_0 = floor(ff_00/(86400/L));
neg1 = ff_0;

gg_0 = floor(gg_00/(86400/L));
neg2 = ff_0;

hh_0 = floor(hh_00/(86400/L));
neg3 = hh_0;

ii_0 = floor(ii_00/(86400/L));
neg4 = hh_0;



% 구간의 시작점과 끝점을 변수에 저장
start_points1 = ff_0(:, 1);
end_points1 = ff_0(:, 2);

start_points2 = ff_0(:, 1);
end_points2 = ff_0(:, 2);

start_points3 = hh_0(:, 1);
end_points3 = hh_0(:, 2);

start_points4 = hh_0(:, 1);
end_points4 = hh_0(:, 2);


% 모든 구간을 정수로 표시해서 작업 공간에 저장
result_kh1 = [];
result_kh2 = [];
result_kh3 = [];
result_kh4 = [];
result_kh5 = [];
result_kh6 = [];
result_kh7 = [];
result_kh8 = [];
result_kh9 = [];
result_kh10 = [];
result_kh11 = [];

result_eh1 = [];
result_eh2 = [];
result_eh3 = [];
result_eh4 = [];
result_eh5 = [];
result_eh6 = [];
result_eh7 = [];
result_eh8 = [];
result_eh9 = [];
result_eh10 = [];
result_eh11 = [];


%% K_Satellite Setting(L0 / H2)

for n = 0:719
integer_intervals_kh1 = [];
integer_intervals_kh2 = [];

  
   for i = 1:length(ff_0)
   % 고궤도 satellite 1 위치 생성
   interval_kh1 = start_points1(i) :end_points1(i) ;
       interval_kh1(interval_kh1 > 719) = interval_kh1(interval_kh1 > 719) - 720; 
       integer_intervals_kh1 = [integer_intervals_kh1, interval_kh1];

   % 고궤도 satellite 2 랜덤 위치 생성
       interval_kh2 = interval_kh1 + n; 
       interval_kh2(interval_kh2 > 719) = interval_kh2(interval_kh2 > 719) - 720;
       integer_intervals_kh2 = [integer_intervals_kh2, interval_kh2];
   end


   sum_intervals_kh = ...
       [integer_intervals_kh1, ...
       integer_intervals_kh2, ...
       ];


unique_intervals_kh = unique(sort(sum_intervals_kh));

difference_kh = length(sum_intervals_kh) - length(unique_intervals_kh) ;

result_kh1 = [result_kh1; n, difference_kh];

end

Cov_kh1 = length(integer_intervals_kh1)/720*100;

%% E_Satellite Setting(L0 / H2)

for n = 0:719
integer_intervals_eh1 = [];
integer_intervals_eh2 = [];

  
   for i = 1:length(hh_0)
   % 고궤도 satellite 1 위치 생성
       interval_eh1 = start_points3(i) :end_points3(i) ;
       interval_eh1(interval_eh1 > 719) = interval_eh1(interval_eh1 > 719) - 720; 
       integer_intervals_eh1 = [integer_intervals_eh1, interval_eh1];

   % 고궤도 satellite 2 랜덤 위치 생성
       interval_eh2 = interval_eh1 + n; 
       interval_eh2(interval_eh2 > 719) = interval_eh2(interval_eh2 > 719) - 720;
       integer_intervals_eh2 = [integer_intervals_eh2, interval_eh2];
   end


   sum_intervals_eh = ...
       [integer_intervals_eh1, ...
       integer_intervals_eh2, ...
       ];


unique_intervals_eh = unique(sort(sum_intervals_eh));

difference_eh = length(sum_intervals_eh) - length(unique_intervals_eh) ;

result_eh1 = [result_eh1; n, difference_eh];

end

Cov_eh1 = length(integer_intervals_eh1)/720*100;

%% 공통

result_sum1 = [result_kh1(:,1), result_kh1(:,2) + result_eh1(:,2)];


rows_with_min_sum1 = result_sum1(result_sum1(:, 2) == min(result_sum1(:, 2)), :);

k_sum1 = rows_with_min_sum1(:,1);

random_index_1 = randi([1, length(k_sum1)]);

random_number_1 = k_sum1(random_index_1);

%% K_renewal 2
integer_intervals_kh1 = [];
integer_intervals_kh2 = [];

   for i = 1:length(ff_0)
   % 고궤도 satellite 1 위치 생성
   interval_kh1 = start_points1(i) :end_points1(i) ;
       interval_kh1(interval_kh1 > 719) = interval_kh1(interval_kh1 > 719) - 720; 
       integer_intervals_kh1 = [integer_intervals_kh1, interval_kh1];

   % 고궤도 satellite 2 위치 재생성
       interval_kh2 = interval_kh1 + random_number_1; 
       interval_kh2(interval_kh2 > 719) = interval_kh2(interval_kh2 > 719) - 720;
       integer_intervals_kh2 = [integer_intervals_kh2, interval_kh2];
   end


   sum_intervals_kh = ...
       [integer_intervals_kh1, ...
       integer_intervals_kh2, ...
       ];


unique_intervals_kh = unique(sort(sum_intervals_kh));

Cov_kh2 = length(unique_intervals_kh) / 720 * 100;

%% E_renewal 2
integer_intervals_eh1 = [];
integer_intervals_eh2 = [];

   for i = 1:length(hh_0)
   % 고궤도 satellite 1 위치 생성
   interval_eh1 = start_points3(i) :end_points3(i) ;
       interval_eh1(interval_eh1 > 719) = interval_eh1(interval_eh1 > 719) - 720; 
       integer_intervals_eh1 = [integer_intervals_eh1, interval_eh1];

   % 고궤도 satellite 2 위치 재생성
       interval_eh2 = interval_eh1 + random_number_1; 
       interval_eh2(interval_eh2 > 719) = interval_eh2(interval_eh2 > 719) - 720;
       integer_intervals_eh2 = [integer_intervals_eh2, interval_eh2];
   end


   sum_intervals_eh = ...
       [integer_intervals_eh1, ...
       integer_intervals_eh2, ...
       ];


unique_intervals_eh = unique(sort(sum_intervals_eh));

Cov_eh2 = length(unique_intervals_eh) / 720 * 100;

%%
%% K_Satellite Setting(L0 / H3)

for n = 0:719
integer_intervals_kh1 = [];
integer_intervals_kh2 = [];
integer_intervals_kh3 = [];

  
   for i = 1:length(ff_0)
   % 고궤도 satellite 1 위치 생성
       interval_kh1 = start_points1(i) :end_points1(i) ;
       interval_kh1(interval_kh1 > 719) = interval_kh1(interval_kh1 > 719) - 720; 
       integer_intervals_kh1 = [integer_intervals_kh1, interval_kh1];

   % 고궤도 satellite 2 위치 생성
       interval_kh2 = interval_kh1 + random_number_1; 
       interval_kh2(interval_kh2 > 719) = interval_kh2(interval_kh2 > 719) - 720;
       integer_intervals_kh2 = [integer_intervals_kh2, interval_kh2];

   % 고궤도 satellite 3 랜덤 위치 생성
       interval_kh3 = interval_kh1 + n; 
       interval_kh3(interval_kh3 > 719) = interval_kh3(interval_kh3 > 719) - 720;
       integer_intervals_kh3 = [integer_intervals_kh3, interval_kh3];
   end


   sum_intervals_kh = ...
       [integer_intervals_kh1, ...
       integer_intervals_kh2, ...
       integer_intervals_kh3, ...
       ];


unique_intervals_kh = unique(sort(sum_intervals_kh));

difference_kh = length(sum_intervals_kh) - length(unique_intervals_kh) ;

result_kh2 = [result_kh2; n, difference_kh];

end

%% E_Satellite Setting(L0 / H3)

for n = 0:719
integer_intervals_eh1 = [];
integer_intervals_eh2 = [];
integer_intervals_eh3 = [];

  
   for i = 1:length(hh_0)
   % 고궤도 satellite 1 위치 생성
       interval_eh1 = start_points3(i) :end_points3(i) ;
       interval_eh1(interval_eh1 > 719) = interval_eh1(interval_eh1 > 719) - 720; 
       integer_intervals_eh1 = [integer_intervals_eh1, interval_eh1];

   % 고궤도 satellite 2 위치 생성
       interval_eh2 = interval_eh1 + random_number_1; 
       interval_eh2(interval_eh2 > 719) = interval_eh2(interval_eh2 > 719) - 720;
       integer_intervals_eh2 = [integer_intervals_eh2, interval_eh2];

   % 고궤도 satellite 3 랜덤 위치 생성
       interval_eh3 = interval_eh1 + n; 
       interval_eh3(interval_eh3 > 719) = interval_eh3(interval_eh3 > 719) - 720;
       integer_intervals_eh3 = [integer_intervals_eh3, interval_eh3];
   end


   sum_intervals_eh = ...
       [integer_intervals_eh1, ...
       integer_intervals_eh2, ...
       integer_intervals_eh3, ...
       ];


unique_intervals_eh = unique(sort(sum_intervals_eh));

difference_eh = length(sum_intervals_eh) - length(unique_intervals_eh) ;

result_eh2 = [result_eh2; n, difference_eh];

end

%% 공통

result_sum2 = [result_kh2(:,1), result_kh2(:,2) + result_eh2(:,2)];


rows_with_min_sum2 = result_sum2(result_sum2(:, 2) == min(result_sum2(:, 2)), :);

k_sum2 = rows_with_min_sum2(:,1);

random_index_2 = randi([1, length(k_sum2)]);

random_number_2 = k_sum2(random_index_2);

%% K_renewal 3
integer_intervals_kh1 = [];
integer_intervals_kh2 = [];
integer_intervals_kh3 = [];

   for i = 1:length(ff_0)
   % 고궤도 satellite 1 위치 생성
   interval_kh1 = start_points1(i) :end_points1(i) ;
       interval_kh1(interval_kh1 > 719) = interval_kh1(interval_kh1 > 719) - 720; 
       integer_intervals_kh1 = [integer_intervals_kh1, interval_kh1];

   % 고궤도 satellite 2 위치 재생성
       interval_kh2 = interval_kh1 + random_number_1; 
       interval_kh2(interval_kh2 > 719) = interval_kh2(interval_kh2 > 719) - 720;
       integer_intervals_kh2 = [integer_intervals_kh2, interval_kh2];

   % 고궤도 satellite 3 위치 재생성
       interval_kh3 = interval_kh1 + random_number_2; 
       interval_kh3(interval_kh3 > 719) = interval_kh3(interval_kh3 > 719) - 720;
       integer_intervals_kh3 = [integer_intervals_kh3, interval_kh3];
   end


   sum_intervals_kh = ...
       [integer_intervals_kh1, ...
       integer_intervals_kh2, ...
       integer_intervals_kh3, ...
       ];


unique_intervals_kh = unique(sort(sum_intervals_kh));

Cov_kh3 = length(unique_intervals_kh) / 720 * 100;

%% E_renewal 3
integer_intervals_eh1 = [];
integer_intervals_eh2 = [];
integer_intervals_eh3 = [];

   for i = 1:length(hh_0)
   % 고궤도 satellite 1 위치 생성
   interval_eh1 = start_points3(i) :end_points3(i) ;
       interval_eh1(interval_eh1 > 719) = interval_eh1(interval_eh1 > 719) - 720; 
       integer_intervals_eh1 = [integer_intervals_eh1, interval_eh1];

   % 고궤도 satellite 2 위치 재생성
       interval_eh2 = interval_eh1 + random_number_1; 
       interval_eh2(interval_eh2 > 719) = interval_eh2(interval_eh2 > 719) - 720;
       integer_intervals_eh2 = [integer_intervals_eh2, interval_eh2];

   % 고궤도 satellite 3 위치 재생성
       interval_eh3 = interval_eh1 + random_number_2; 
       interval_eh3(interval_eh3 > 719) = interval_eh3(interval_eh3 > 719) - 720;
       integer_intervals_eh3 = [integer_intervals_eh3, interval_eh3];
   end


   sum_intervals_eh = ...
       [integer_intervals_eh1, ...
       integer_intervals_eh2, ...
       integer_intervals_eh3, ...
       ];


unique_intervals_eh = unique(sort(sum_intervals_eh));

Cov_eh3 = length(unique_intervals_eh) / 720 * 100;

%%
%% K_Satellite Setting(L0 / H4)

for n = 0:719
integer_intervals_kh1 = [];
integer_intervals_kh2 = [];
integer_intervals_kh3 = [];
integer_intervals_kh4 = [];

  
   for i = 1:length(ff_0)
   % 고궤도 satellite 1 위치 생성
       interval_kh1 = start_points1(i) :end_points1(i) ;
       interval_kh1(interval_kh1 > 719) = interval_kh1(interval_kh1 > 719) - 720; 
       integer_intervals_kh1 = [integer_intervals_kh1, interval_kh1];

   % 고궤도 satellite 2 위치 생성
       interval_kh2 = interval_kh1 + random_number_1; 
       interval_kh2(interval_kh2 > 719) = interval_kh2(interval_kh2 > 719) - 720;
       integer_intervals_kh2 = [integer_intervals_kh2, interval_kh2];

   % 고궤도 satellite 3 위치 생성
       interval_kh3 = interval_kh1 + random_number_2; 
       interval_kh3(interval_kh3 > 719) = interval_kh3(interval_kh3 > 719) - 720;
       integer_intervals_kh3 = [integer_intervals_kh3, interval_kh3];

   % 고궤도 satellite 4 랜덤 위치 생성
       interval_kh4 = interval_kh1 + n; 
       interval_kh4(interval_kh4 > 719) = interval_kh4(interval_kh4 > 719) - 720;
       integer_intervals_kh4 = [integer_intervals_kh4, interval_kh4];
   end


   sum_intervals_kh = ...
       [integer_intervals_kh1, ...
       integer_intervals_kh2, ...
       integer_intervals_kh3, ...
       integer_intervals_kh4, ...
       ];


unique_intervals_kh = unique(sort(sum_intervals_kh));

difference_kh = length(sum_intervals_kh) - length(unique_intervals_kh) ;

result_kh3 = [result_kh3; n, difference_kh];

end

%% E_Satellite Setting(L0 / H4)

for n = 0:719
integer_intervals_eh1 = [];
integer_intervals_eh2 = [];
integer_intervals_eh3 = [];
integer_intervals_eh4 = [];

  
   for i = 1:length(hh_0)
   % 고궤도 satellite 1 위치 생성
       interval_eh1 = start_points3(i) :end_points3(i) ;
       interval_eh1(interval_eh1 > 719) = interval_eh1(interval_eh1 > 719) - 720; 
       integer_intervals_eh1 = [integer_intervals_eh1, interval_eh1];

   % 고궤도 satellite 2 위치 생성
       interval_eh2 = interval_eh1 + random_number_1; 
       interval_eh2(interval_eh2 > 719) = interval_eh2(interval_eh2 > 719) - 720;
       integer_intervals_eh2 = [integer_intervals_eh2, interval_eh2];

   % 고궤도 satellite 3 위치 생성
       interval_eh3 = interval_eh1 + random_number_2; 
       interval_eh3(interval_eh3 > 719) = interval_eh3(interval_eh3 > 719) - 720;
       integer_intervals_eh3 = [integer_intervals_eh3, interval_eh3];

   % 고궤도 satellite 4 랜덤 위치 생성
       interval_eh4 = interval_eh1 + n; 
       interval_eh4(interval_eh4 > 719) = interval_eh4(interval_eh4 > 719) - 720;
       integer_intervals_eh4 = [integer_intervals_eh4, interval_eh4];
   end


   sum_intervals_eh = ...
       [integer_intervals_eh1, ...
       integer_intervals_eh2, ...
       integer_intervals_eh3, ...
       integer_intervals_eh4, ...
       ];


unique_intervals_eh = unique(sort(sum_intervals_eh));

difference_eh = length(sum_intervals_eh) - length(unique_intervals_eh) ;

result_eh3 = [result_eh3; n, difference_eh];

end

%% 공통

result_sum3 = [result_kh3(:,1), result_kh3(:,2) + result_eh3(:,2)];


rows_with_min_sum3 = result_sum3(result_sum3(:, 2) == min(result_sum3(:, 2)), :);

k_sum3 = rows_with_min_sum3(:,1);

random_index_3 = randi([1, length(k_sum3)]);

random_number_3 = k_sum3(random_index_3);

%% K_renewal 4
integer_intervals_kh1 = [];
integer_intervals_kh2 = [];
integer_intervals_kh3 = [];
integer_intervals_kh4 = [];

   for i = 1:length(ff_0)
   % 고궤도 satellite 1 위치 생성
       interval_kh1 = start_points1(i) :end_points1(i) ;
       interval_kh1(interval_kh1 > 719) = interval_kh1(interval_kh1 > 719) - 720; 
       integer_intervals_kh1 = [integer_intervals_kh1, interval_kh1];

   % 고궤도 satellite 2 위치 재생성
       interval_kh2 = interval_kh1 + random_number_1; 
       interval_kh2(interval_kh2 > 719) = interval_kh2(interval_kh2 > 719) - 720;
       integer_intervals_kh2 = [integer_intervals_kh2, interval_kh2];

   % 고궤도 satellite 3 위치 재생성
       interval_kh3 = interval_kh1 + random_number_2; 
       interval_kh3(interval_kh3 > 719) = interval_kh3(interval_kh3 > 719) - 720;
       integer_intervals_kh3 = [integer_intervals_kh3, interval_kh3];

   % 고궤도 satellite 4 위치 재생성
       interval_kh4 = interval_kh1 + random_number_3; 
       interval_kh4(interval_kh4 > 719) = interval_kh4(interval_kh4 > 719) - 720;
       integer_intervals_kh4 = [integer_intervals_kh4, interval_kh4];
   end


   sum_intervals_kh = ...
       [integer_intervals_kh1, ...
       integer_intervals_kh2, ...
       integer_intervals_kh3, ...
       integer_intervals_kh4, ...
       ];


unique_intervals_kh = unique(sort(sum_intervals_kh));

Cov_kh4 = length(unique_intervals_kh) / 720 * 100;

%% E_renewal 4
integer_intervals_eh1 = [];
integer_intervals_eh2 = [];
integer_intervals_eh3 = [];
integer_intervals_eh4 = [];

   for i = 1:length(hh_0)
   % 고궤도 satellite 1 위치 생성
   interval_eh1 = start_points3(i) :end_points3(i) ;
       interval_eh1(interval_eh1 > 719) = interval_eh1(interval_eh1 > 719) - 720; 
       integer_intervals_eh1 = [integer_intervals_eh1, interval_eh1];

   % 고궤도 satellite 2 위치 재생성
       interval_eh2 = interval_eh1 + random_number_1; 
       interval_eh2(interval_eh2 > 719) = interval_eh2(interval_eh2 > 719) - 720;
       integer_intervals_eh2 = [integer_intervals_eh2, interval_eh2];

   % 고궤도 satellite 3 위치 재생성
       interval_eh3 = interval_eh1 + random_number_2; 
       interval_eh3(interval_eh3 > 719) = interval_eh3(interval_eh3 > 719) - 720;
       integer_intervals_eh3 = [integer_intervals_eh3, interval_eh3];

   % 고궤도 satellite 4 위치 재생성
       interval_eh4 = interval_eh1 + random_number_3; 
       interval_eh4(interval_eh4 > 719) = interval_eh4(interval_eh4 > 719) - 720;
       integer_intervals_eh4 = [integer_intervals_eh4, interval_eh4];
   end


   sum_intervals_eh = ...
       [integer_intervals_eh1, ...
       integer_intervals_eh2, ...
       integer_intervals_eh3, ...
       integer_intervals_eh4, ...
       ];


unique_intervals_eh = unique(sort(sum_intervals_eh));

Cov_eh4 = length(unique_intervals_eh) / 720 * 100;

%%
%% K_Satellite Setting(L0 / H5)

for n = 0:719
integer_intervals_kh1 = [];
integer_intervals_kh2 = [];
integer_intervals_kh3 = [];
integer_intervals_kh4 = [];
integer_intervals_kh5 = [];

  
   for i = 1:length(ff_0)
   % 고궤도 satellite 1 위치 생성
       interval_kh1 = start_points1(i) :end_points1(i) ;
       interval_kh1(interval_kh1 > 719) = interval_kh1(interval_kh1 > 719) - 720; 
       integer_intervals_kh1 = [integer_intervals_kh1, interval_kh1];

   % 고궤도 satellite 2 위치 생성5
       interval_kh2 = interval_kh1 + random_number_1; 
       interval_kh2(interval_kh2 > 719) = interval_kh2(interval_kh2 > 719) - 720;
       integer_intervals_kh2 = [integer_intervals_kh2, interval_kh2];

   % 고궤도 satellite 3 위치 생성
       interval_kh3 = interval_kh1 + random_number_2; 
       interval_kh3(interval_kh3 > 719) = interval_kh3(interval_kh3 > 719) - 720;
       integer_intervals_kh3 = [integer_intervals_kh3, interval_kh3];

   % 고궤도 satellite 4 위치 생성
       interval_kh4 = interval_kh1 + random_number_3; 
       interval_kh4(interval_kh4 > 719) = interval_kh4(interval_kh4 > 719) - 720;
       integer_intervals_kh4 = [integer_intervals_kh4, interval_kh4];

   % 고궤도 satellite 5 랜덤 위치 생성
       interval_kh5 = interval_kh1 + n; 
       interval_kh5(interval_kh5 > 719) = interval_kh5(interval_kh5 > 719) - 720;
       integer_intervals_kh5 = [integer_intervals_kh5, interval_kh5];
   end


   sum_intervals_kh = ...
       [integer_intervals_kh1, ...
       integer_intervals_kh2, ...
       integer_intervals_kh3, ...
       integer_intervals_kh4, ...
       integer_intervals_kh5, ...
       ];


unique_intervals_kh = unique(sort(sum_intervals_kh));

difference_kh = length(sum_intervals_kh) - length(unique_intervals_kh) ;

result_kh4 = [result_kh4; n, difference_kh];

end

%% E_Satellite Setting(L0 / H5)

for n = 0:719
integer_intervals_eh1 = [];
integer_intervals_eh2 = [];
integer_intervals_eh3 = [];
integer_intervals_eh4 = [];
integer_intervals_eh5 = [];

  
   for i = 1:length(hh_0)
   % 고궤도 satellite 1 위치 생성
       interval_eh1 = start_points3(i) :end_points3(i) ;
       interval_eh1(interval_eh1 > 719) = interval_eh1(interval_eh1 > 719) - 720; 
       integer_intervals_eh1 = [integer_intervals_eh1, interval_eh1];

   % 고궤도 satellite 2 위치 생성
       interval_eh2 = interval_eh1 + random_number_1; 
       interval_eh2(interval_eh2 > 719) = interval_eh2(interval_eh2 > 719) - 720;
       integer_intervals_eh2 = [integer_intervals_eh2, interval_eh2];

   % 고궤도 satellite 3 위치 생성
       interval_eh3 = interval_eh1 + random_number_2; 
       interval_eh3(interval_eh3 > 719) = interval_eh3(interval_eh3 > 719) - 720;
       integer_intervals_eh3 = [integer_intervals_eh3, interval_eh3];

   % 고궤도 satellite 4 위치 생성
       interval_eh4 = interval_eh1 + random_number_3; 
       interval_eh4(interval_eh4 > 719) = interval_eh4(interval_eh4 > 719) - 720;
       integer_intervals_eh4 = [integer_intervals_eh4, interval_eh4];

   % 고궤도 satellite 5 랜덤 위치 생성
       interval_eh5 = interval_eh1 + n; 
       interval_eh5(interval_eh5 > 719) = interval_eh5(interval_eh5 > 719) - 720;
       integer_intervals_eh5 = [integer_intervals_eh5, interval_eh5];
   end


   sum_intervals_eh = ...
       [integer_intervals_eh1, ...
       integer_intervals_eh2, ...
       integer_intervals_eh3, ...
       integer_intervals_eh4, ...
       integer_intervals_eh5, ...
       ];


unique_intervals_eh = unique(sort(sum_intervals_eh));

difference_eh = length(sum_intervals_eh) - length(unique_intervals_eh) ;

result_eh4 = [result_eh4; n, difference_eh];

end

%% 공통

result_sum4 = [result_kh4(:,1), result_kh4(:,2) + result_eh4(:,2)];


rows_with_min_sum4 = result_sum4(result_sum4(:, 2) == min(result_sum4(:, 2)), :);

k_sum4 = rows_with_min_sum4(:,1);

random_index_4 = randi([1, length(k_sum4)]);

random_number_4 = k_sum4(random_index_4);

%% K_renewal 5
integer_intervals_kh1 = [];
integer_intervals_kh2 = [];
integer_intervals_kh3 = [];
integer_intervals_kh4 = [];
integer_intervals_kh5 = [];

   for i = 1:length(ff_0)
   % 고궤도 satellite 1 위치 생성
       interval_kh1 = start_points1(i) :end_points1(i) ;
       interval_kh1(interval_kh1 > 719) = interval_kh1(interval_kh1 > 719) - 7320; 
       integer_intervals_kh1 = [integer_intervals_kh1, interval_kh1];

   % 고궤도 satellite 2 위치 재생성
       interval_kh2 = interval_kh1 + random_number_1; 
       interval_kh2(interval_kh2 > 719) = interval_kh2(interval_kh2 > 719) - 720;
       integer_intervals_kh2 = [integer_intervals_kh2, interval_kh2];

   % 고궤도 satellite 3 위치 재생성
       interval_kh3 = interval_kh1 + random_number_2; 
       interval_kh3(interval_kh3 > 719) = interval_kh3(interval_kh3 > 719) - 720;
       integer_intervals_kh3 = [integer_intervals_kh3, interval_kh3];

   % 고궤도 satellite 4 위치 재생성
       interval_kh4 = interval_kh1 + random_number_3; 
       interval_kh4(interval_kh4 > 719) = interval_kh4(interval_kh4 > 719) - 720;
       integer_intervals_kh4 = [integer_intervals_kh4, interval_kh4];

   % 고궤도 satellite 5 위치 재생성
       interval_kh5 = interval_kh1 + random_number_4; 
       interval_kh5(interval_kh5 > 719) = interval_kh5(interval_kh5 > 719) - 720;
       integer_intervals_kh5 = [integer_intervals_kh5, interval_kh5];
   end


   sum_intervals_kh = ...
       [integer_intervals_kh1, ...
       integer_intervals_kh2, ...
       integer_intervals_kh3, ...
       integer_intervals_kh4, ...
       integer_intervals_kh5, ...
       ];


unique_intervals_kh = unique(sort(sum_intervals_kh));

Cov_kh5 = length(unique_intervals_kh) / 720 * 100;

%% E_renewal 5
integer_intervals_eh1 = [];
integer_intervals_eh2 = [];
integer_intervals_eh3 = [];
integer_intervals_eh4 = [];
integer_intervals_eh5 = [];

   for i = 1:length(hh_0)
   % 고궤도 satellite 1 위치 생성
   interval_eh1 = start_points3(i) :end_points3(i) ;
       interval_eh1(interval_eh1 > 719) = interval_eh1(interval_eh1 > 719) - 720; 
       integer_intervals_eh1 = [integer_intervals_eh1, interval_eh1];

   % 고궤도 satellite 2 위치 재생성
       interval_eh2 = interval_eh1 + random_number_1; 
       interval_eh2(interval_eh2 > 719) = interval_eh2(interval_eh2 > 719) - 720;
       integer_intervals_eh2 = [integer_intervals_eh2, interval_eh2];

   % 고궤도 satellite 3 위치 재생성
       interval_eh3 = interval_eh1 + random_number_2; 
       interval_eh3(interval_eh3 > 719) = interval_eh3(interval_eh3 > 719) - 720;
       integer_intervals_eh3 = [integer_intervals_eh3, interval_eh3];

   % 고궤도 satellite 4 위치 재생성
       interval_eh4 = interval_eh1 + random_number_3; 
       interval_eh4(interval_eh4 > 719) = interval_eh4(interval_eh4 > 719) - 720;
       integer_intervals_eh4 = [integer_intervals_eh4, interval_eh4];

   % 고궤도 satellite 5 위치 재생성
       interval_eh5 = interval_eh1 + random_number_4; 
       interval_eh5(interval_eh5 > 719) = interval_eh5(interval_eh5 > 719) - 720;
       integer_intervals_eh5 = [integer_intervals_eh5, interval_eh5];
   end


   sum_intervals_eh = ...
       [integer_intervals_eh1, ...
       integer_intervals_eh2, ...
       integer_intervals_eh3, ...
       integer_intervals_eh4, ...
       integer_intervals_eh5, ...
       ];


unique_intervals_eh = unique(sort(sum_intervals_eh));

Cov_eh5 = length(unique_intervals_eh) / 720 * 100;


%%
%% K_Satellite Setting(L0 / H6)

for n = 0:719
integer_intervals_kh1 = [];
integer_intervals_kh2 = [];
integer_intervals_kh3 = [];
integer_intervals_kh4 = [];
integer_intervals_kh5 = [];
integer_intervals_kh6 = [];

  
   for i = 1:length(ff_0)
   % 고궤도 satellite 1 위치 생성
       interval_kh1 = start_points1(i) :end_points1(i) ;
       interval_kh1(interval_kh1 > 719) = interval_kh1(interval_kh1 > 719) - 720; 
       integer_intervals_kh1 = [integer_intervals_kh1, interval_kh1];

   % 고궤도 satellite 2 위치 생성5
       interval_kh2 = interval_kh1 + random_number_1; 
       interval_kh2(interval_kh2 > 719) = interval_kh2(interval_kh2 > 719) - 720;
       integer_intervals_kh2 = [integer_intervals_kh2, interval_kh2];

   % 고궤도 satellite 3 위치 생성
       interval_kh3 = interval_kh1 + random_number_2; 
       interval_kh3(interval_kh3 > 719) = interval_kh3(interval_kh3 > 719) - 720;
       integer_intervals_kh3 = [integer_intervals_kh3, interval_kh3];

   % 고궤도 satellite 4 위치 생성
       interval_kh4 = interval_kh1 + random_number_3; 
       interval_kh4(interval_kh4 > 719) = interval_kh4(interval_kh4 > 719) - 720;
       integer_intervals_kh4 = [integer_intervals_kh4, interval_kh4];

   % 고궤도 satellite 5 위치 생성
       interval_kh5 = interval_kh1 + random_number_4; 
       interval_kh5(interval_kh5 > 719) = interval_kh5(interval_kh5 > 719) - 720;
       integer_intervals_kh5 = [integer_intervals_kh5, interval_kh5];

   % 고궤도 satellite 6 랜덤 위치 생성
       interval_kh6 = interval_kh1 + n; 
       interval_kh6(interval_kh6 > 719) = interval_kh6(interval_kh6 > 719) - 720;
       integer_intervals_kh6 = [integer_intervals_kh6, interval_kh6];
   end


   sum_intervals_kh = ...
       [integer_intervals_kh1, ...
       integer_intervals_kh2, ...
       integer_intervals_kh3, ...
       integer_intervals_kh4, ...
       integer_intervals_kh5, ...
       integer_intervals_kh6, ...
       ];


unique_intervals_kh = unique(sort(sum_intervals_kh));

difference_kh = length(sum_intervals_kh) - length(unique_intervals_kh) ;

result_kh5 = [result_kh5; n, difference_kh];

end

%% E_Satellite Setting(L0 / H6)

for n = 0:719
integer_intervals_eh1 = [];
integer_intervals_eh2 = [];
integer_intervals_eh3 = [];
integer_intervals_eh4 = [];
integer_intervals_eh5 = [];
integer_intervals_eh6 = [];

  
   for i = 1:length(hh_0)
   % 고궤도 satellite 1 위치 생성
       interval_eh1 = start_points3(i) :end_points3(i) ;
       interval_eh1(interval_eh1 > 719) = interval_eh1(interval_eh1 > 719) - 720; 
       integer_intervals_eh1 = [integer_intervals_eh1, interval_eh1];

   % 고궤도 satellite 2 위치 생성
       interval_eh2 = interval_eh1 + random_number_1; 
       interval_eh2(interval_eh2 > 719) = interval_eh2(interval_eh2 > 719) - 720;
       integer_intervals_eh2 = [integer_intervals_eh2, interval_eh2];

   % 고궤도 satellite 3 위치 생성
       interval_eh3 = interval_eh1 + random_number_2; 
       interval_eh3(interval_eh3 > 719) = interval_eh3(interval_eh3 > 719) - 720;
       integer_intervals_eh3 = [integer_intervals_eh3, interval_eh3];

   % 고궤도 satellite 4 위치 생성
       interval_eh4 = interval_eh1 + random_number_3; 
       interval_eh4(interval_eh4 > 719) = interval_eh4(interval_eh4 > 719) - 720;
       integer_intervals_eh4 = [integer_intervals_eh4, interval_eh4];

   % 고궤도 satellite 5 위치 생성
       interval_eh5 = interval_eh1 + random_number_4; 
       interval_eh5(interval_eh5 > 719) = interval_eh5(interval_eh5 > 719) - 720;
       integer_intervals_eh5 = [integer_intervals_eh5, interval_eh5];

   % 고궤도 satellite 6 랜덤 위치 생성
       interval_eh6 = interval_eh1 + n; 
       interval_eh6(interval_eh6 > 719) = interval_eh6(interval_eh6 > 719) - 720;
       integer_intervals_eh6 = [integer_intervals_eh6, interval_eh6];
   end


   sum_intervals_eh = ...
       [integer_intervals_eh1, ...
       integer_intervals_eh2, ...
       integer_intervals_eh3, ...
       integer_intervals_eh4, ...
       integer_intervals_eh5, ...
       integer_intervals_eh6, ...
       ];


unique_intervals_eh = unique(sort(sum_intervals_eh));

difference_eh = length(sum_intervals_eh) - length(unique_intervals_eh) ;

result_eh5 = [result_eh5; n, difference_eh];

end

%% 공통

result_sum5 = [result_kh5(:,1), result_kh5(:,2) + result_eh5(:,2)];


rows_with_min_sum5 = result_sum5(result_sum5(:, 2) == min(result_sum5(:, 2)), :);

k_sum5 = rows_with_min_sum5(:,1);

random_index_5 = randi([1, length(k_sum5)]);

random_number_5 = k_sum5(random_index_5);

%% K_renewal 6
integer_intervals_kh1 = [];
integer_intervals_kh2 = [];
integer_intervals_kh3 = [];
integer_intervals_kh4 = [];
integer_intervals_kh5 = [];
integer_intervals_kh6 = [];

   for i = 1:length(ff_0)
   % 고궤도 satellite 1 위치 생성
       interval_kh1 = start_points1(i) :end_points1(i) ;
       interval_kh1(interval_kh1 > 719) = interval_kh1(interval_kh1 > 719) - 7320; 
       integer_intervals_kh1 = [integer_intervals_kh1, interval_kh1];

   % 고궤도 satellite 2 위치 재생성
       interval_kh2 = interval_kh1 + random_number_1; 
       interval_kh2(interval_kh2 > 719) = interval_kh2(interval_kh2 > 719) - 720;
       integer_intervals_kh2 = [integer_intervals_kh2, interval_kh2];

   % 고궤도 satellite 3 위치 재생성
       interval_kh3 = interval_kh1 + random_number_2; 
       interval_kh3(interval_kh3 > 719) = interval_kh3(interval_kh3 > 719) - 720;
       integer_intervals_kh3 = [integer_intervals_kh3, interval_kh3];

   % 고궤도 satellite 4 위치 재생성
       interval_kh4 = interval_kh1 + random_number_3; 
       interval_kh4(interval_kh4 > 719) = interval_kh4(interval_kh4 > 719) - 720;
       integer_intervals_kh4 = [integer_intervals_kh4, interval_kh4];

   % 고궤도 satellite 5 위치 재생성
       interval_kh5 = interval_kh1 + random_number_4; 
       interval_kh5(interval_kh5 > 719) = interval_kh5(interval_kh5 > 719) - 720;
       integer_intervals_kh5 = [integer_intervals_kh5, interval_kh5];

   % 고궤도 satellite 6 위치 재생성
       interval_kh6 = interval_kh1 + random_number_5; 
       interval_kh6(interval_kh6 > 719) = interval_kh6(interval_kh6 > 719) - 720;
       integer_intervals_kh6 = [integer_intervals_kh6, interval_kh6];
   end


   sum_intervals_kh = ...
       [integer_intervals_kh1, ...
       integer_intervals_kh2, ...
       integer_intervals_kh3, ...
       integer_intervals_kh4, ...
       integer_intervals_kh5, ...
       integer_intervals_kh6, ...
       ];


unique_intervals_kh = unique(sort(sum_intervals_kh));

Cov_kh6 = length(unique_intervals_kh) / 720 * 100;

%% E_renewal 6
integer_intervals_eh1 = [];
integer_intervals_eh2 = [];
integer_intervals_eh3 = [];
integer_intervals_eh4 = [];
integer_intervals_eh5 = [];
integer_intervals_eh6 = [];

   for i = 1:length(hh_0)
   % 고궤도 satellite 1 위치 생성
   interval_eh1 = start_points3(i) :end_points3(i) ;
       interval_eh1(interval_eh1 > 719) = interval_eh1(interval_eh1 > 719) - 720; 
       integer_intervals_eh1 = [integer_intervals_eh1, interval_eh1];

   % 고궤도 satellite 2 위치 재생성
       interval_eh2 = interval_eh1 + random_number_1; 
       interval_eh2(interval_eh2 > 719) = interval_eh2(interval_eh2 > 719) - 720;
       integer_intervals_eh2 = [integer_intervals_eh2, interval_eh2];

   % 고궤도 satellite 3 위치 재생성
       interval_eh3 = interval_eh1 + random_number_2; 
       interval_eh3(interval_eh3 > 719) = interval_eh3(interval_eh3 > 719) - 720;
       integer_intervals_eh3 = [integer_intervals_eh3, interval_eh3];

   % 고궤도 satellite 4 위치 재생성
       interval_eh4 = interval_eh1 + random_number_3; 
       interval_eh4(interval_eh4 > 719) = interval_eh4(interval_eh4 > 719) - 720;
       integer_intervals_eh4 = [integer_intervals_eh4, interval_eh4];

   % 고궤도 satellite 5 위치 재생성
       interval_eh5 = interval_eh1 + random_number_4; 
       interval_eh5(interval_eh5 > 719) = interval_eh5(interval_eh5 > 719) - 720;
       integer_intervals_eh5 = [integer_intervals_eh5, interval_eh5];

   % 고궤도 satellite 6 위치 재생성
       interval_eh6 = interval_eh1 + random_number_5; 
       interval_eh6(interval_eh6 > 719) = interval_eh6(interval_eh6 > 719) - 720;
       integer_intervals_eh6 = [integer_intervals_eh6, interval_eh6];
   end


   sum_intervals_eh = ...
       [integer_intervals_eh1, ...
       integer_intervals_eh2, ...
       integer_intervals_eh3, ...
       integer_intervals_eh4, ...
       integer_intervals_eh5, ...
       integer_intervals_eh6, ...
       ];


unique_intervals_eh = unique(sort(sum_intervals_eh));

Cov_eh6 = length(unique_intervals_eh) / 720 * 100;


%%
%% K_Satellite Setting(L0 / H7)

for n = 0:719
integer_intervals_kh1 = [];
integer_intervals_kh2 = [];
integer_intervals_kh3 = [];
integer_intervals_kh4 = [];
integer_intervals_kh5 = [];
integer_intervals_kh6 = [];
integer_intervals_kh7 = [];

  
   for i = 1:length(ff_0)
   % 고궤도 satellite 1 위치 생성
       interval_kh1 = start_points1(i) :end_points1(i) ;
       interval_kh1(interval_kh1 > 719) = interval_kh1(interval_kh1 > 719) - 720; 
       integer_intervals_kh1 = [integer_intervals_kh1, interval_kh1];

   % 고궤도 satellite 2 위치 생성5
       interval_kh2 = interval_kh1 + random_number_1; 
       interval_kh2(interval_kh2 > 719) = interval_kh2(interval_kh2 > 719) - 720;
       integer_intervals_kh2 = [integer_intervals_kh2, interval_kh2];

   % 고궤도 satellite 3 위치 생성
       interval_kh3 = interval_kh1 + random_number_2; 
       interval_kh3(interval_kh3 > 719) = interval_kh3(interval_kh3 > 719) - 720;
       integer_intervals_kh3 = [integer_intervals_kh3, interval_kh3];

   % 고궤도 satellite 4 위치 생성
       interval_kh4 = interval_kh1 + random_number_3; 
       interval_kh4(interval_kh4 > 719) = interval_kh4(interval_kh4 > 719) - 720;
       integer_intervals_kh4 = [integer_intervals_kh4, interval_kh4];

   % 고궤도 satellite 5 위치 생성
       interval_kh5 = interval_kh1 + random_number_4; 
       interval_kh5(interval_kh5 > 719) = interval_kh5(interval_kh5 > 719) - 720;
       integer_intervals_kh5 = [integer_intervals_kh5, interval_kh5];

   % 고궤도 satellite 6 위치 생성
       interval_kh6 = interval_kh1 + random_number_5; 
       interval_kh6(interval_kh6 > 719) = interval_kh6(interval_kh6 > 719) - 720;
       integer_intervals_kh6 = [integer_intervals_kh6, interval_kh6];

   % 고궤도 satellite 7 랜덤 위치 생성
       interval_kh7 = interval_kh1 + n; 
       interval_kh7(interval_kh7 > 719) = interval_kh7(interval_kh7 > 719) - 720;
       integer_intervals_kh7 = [integer_intervals_kh7, interval_kh7];
   end


   sum_intervals_kh = ...
       [integer_intervals_kh1, ...
       integer_intervals_kh2, ...
       integer_intervals_kh3, ...
       integer_intervals_kh4, ...
       integer_intervals_kh5, ...
       integer_intervals_kh6, ...
       integer_intervals_kh7, ...
       ];


unique_intervals_kh = unique(sort(sum_intervals_kh));

difference_kh = length(sum_intervals_kh) - length(unique_intervals_kh) ;

result_kh6 = [result_kh6; n, difference_kh];

end

%% E_Satellite Setting(L0 / H7)

for n = 0:719
integer_intervals_eh1 = [];
integer_intervals_eh2 = [];
integer_intervals_eh3 = [];
integer_intervals_eh4 = [];
integer_intervals_eh5 = [];
integer_intervals_eh6 = [];
integer_intervals_eh7 = [];

  
   for i = 1:length(hh_0)
   % 고궤도 satellite 1 위치 생성
       interval_eh1 = start_points3(i) :end_points3(i) ;
       interval_eh1(interval_eh1 > 719) = interval_eh1(interval_eh1 > 719) - 720; 
       integer_intervals_eh1 = [integer_intervals_eh1, interval_eh1];

   % 고궤도 satellite 2 위치 생성
       interval_eh2 = interval_eh1 + random_number_1; 
       interval_eh2(interval_eh2 > 719) = interval_eh2(interval_eh2 > 719) - 720;
       integer_intervals_eh2 = [integer_intervals_eh2, interval_eh2];

   % 고궤도 satellite 3 위치 생성
       interval_eh3 = interval_eh1 + random_number_2; 
       interval_eh3(interval_eh3 > 719) = interval_eh3(interval_eh3 > 719) - 720;
       integer_intervals_eh3 = [integer_intervals_eh3, interval_eh3];

   % 고궤도 satellite 4 위치 생성
       interval_eh4 = interval_eh1 + random_number_3; 
       interval_eh4(interval_eh4 > 719) = interval_eh4(interval_eh4 > 719) - 720;
       integer_intervals_eh4 = [integer_intervals_eh4, interval_eh4];

   % 고궤도 satellite 5 위치 생성
       interval_eh5 = interval_eh1 + random_number_4; 
       interval_eh5(interval_eh5 > 719) = interval_eh5(interval_eh5 > 719) - 720;
       integer_intervals_eh5 = [integer_intervals_eh5, interval_eh5];

   % 고궤도 satellite 6 위치 생성
       interval_eh6 = interval_eh1 + random_number_5; 
       interval_eh6(interval_eh6 > 719) = interval_eh6(interval_eh6 > 719) - 720;
       integer_intervals_eh6 = [integer_intervals_eh6, interval_eh6];

   % 고궤도 satellite 7 랜덤 위치 생성
       interval_eh7 = interval_eh1 + n; 
       interval_eh7(interval_eh7 > 719) = interval_eh7(interval_eh7 > 719) - 720;
       integer_intervals_eh7 = [integer_intervals_eh7, interval_eh7];
   end


   sum_intervals_eh = ...
       [integer_intervals_eh1, ...
       integer_intervals_eh2, ...
       integer_intervals_eh3, ...
       integer_intervals_eh4, ...
       integer_intervals_eh5, ...
       integer_intervals_eh6, ...
       integer_intervals_eh7, ...
       ];


unique_intervals_eh = unique(sort(sum_intervals_eh));

difference_eh = length(sum_intervals_eh) - length(unique_intervals_eh) ;

result_eh6 = [result_eh6; n, difference_eh];

end

%% 공통

result_sum6 = [result_kh6(:,1), result_kh6(:,2) + result_eh6(:,2)];


rows_with_min_sum6 = result_sum6(result_sum6(:, 2) == min(result_sum6(:, 2)), :);

k_sum6 = rows_with_min_sum6(:,1);

random_index_6 = randi([1, length(k_sum6)]);

random_number_6 = k_sum6(random_index_6);

%% K_renewal 7
integer_intervals_kh1 = [];
integer_intervals_kh2 = [];
integer_intervals_kh3 = [];
integer_intervals_kh4 = [];
integer_intervals_kh5 = [];
integer_intervals_kh6 = [];
integer_intervals_kh7 = [];

   for i = 1:length(ff_0)
   % 고궤도 satellite 1 위치 생성
       interval_kh1 = start_points1(i) :end_points1(i) ;
       interval_kh1(interval_kh1 > 719) = interval_kh1(interval_kh1 > 719) - 7320; 
       integer_intervals_kh1 = [integer_intervals_kh1, interval_kh1];

   % 고궤도 satellite 2 위치 재생성
       interval_kh2 = interval_kh1 + random_number_1; 
       interval_kh2(interval_kh2 > 719) = interval_kh2(interval_kh2 > 719) - 720;
       integer_intervals_kh2 = [integer_intervals_kh2, interval_kh2];

   % 고궤도 satellite 3 위치 재생성
       interval_kh3 = interval_kh1 + random_number_2; 
       interval_kh3(interval_kh3 > 719) = interval_kh3(interval_kh3 > 719) - 720;
       integer_intervals_kh3 = [integer_intervals_kh3, interval_kh3];

   % 고궤도 satellite 4 위치 재생성
       interval_kh4 = interval_kh1 + random_number_3; 
       interval_kh4(interval_kh4 > 719) = interval_kh4(interval_kh4 > 719) - 720;
       integer_intervals_kh4 = [integer_intervals_kh4, interval_kh4];

   % 고궤도 satellite 5 위치 재생성
       interval_kh5 = interval_kh1 + random_number_4; 
       interval_kh5(interval_kh5 > 719) = interval_kh5(interval_kh5 > 719) - 720;
       integer_intervals_kh5 = [integer_intervals_kh5, interval_kh5];

   % 고궤도 satellite 6 위치 재생성
       interval_kh6 = interval_kh1 + random_number_5; 
       interval_kh6(interval_kh6 > 719) = interval_kh6(interval_kh6 > 719) - 720;
       integer_intervals_kh6 = [integer_intervals_kh6, interval_kh6];

   % 고궤도 satellite 7 위치 재생성
       interval_kh7 = interval_kh1 + random_number_6; 
       interval_kh7(interval_kh7 > 719) = interval_kh7(interval_kh7 > 719) - 720;
       integer_intervals_kh7 = [integer_intervals_kh7, interval_kh7];
   end


   sum_intervals_kh = ...
       [integer_intervals_kh1, ...
       integer_intervals_kh2, ...
       integer_intervals_kh3, ...
       integer_intervals_kh4, ...
       integer_intervals_kh5, ...
       integer_intervals_kh6, ...
       integer_intervals_kh7, ...
       ];


unique_intervals_kh = unique(sort(sum_intervals_kh));

Cov_kh7 = length(unique_intervals_kh) / 720 * 100;

%% E_renewal 7
integer_intervals_eh1 = [];
integer_intervals_eh2 = [];
integer_intervals_eh3 = [];
integer_intervals_eh4 = [];
integer_intervals_eh5 = [];
integer_intervals_eh6 = [];
integer_intervals_eh7 = [];

   for i = 1:length(hh_0)
   % 고궤도 satellite 1 위치 생성
   interval_eh1 = start_points3(i) :end_points3(i) ;
       interval_eh1(interval_eh1 > 719) = interval_eh1(interval_eh1 > 719) - 720; 
       integer_intervals_eh1 = [integer_intervals_eh1, interval_eh1];

   % 고궤도 satellite 2 위치 재생성
       interval_eh2 = interval_eh1 + random_number_1; 
       interval_eh2(interval_eh2 > 719) = interval_eh2(interval_eh2 > 719) - 720;
       integer_intervals_eh2 = [integer_intervals_eh2, interval_eh2];

   % 고궤도 satellite 3 위치 재생성
       interval_eh3 = interval_eh1 + random_number_2; 
       interval_eh3(interval_eh3 > 719) = interval_eh3(interval_eh3 > 719) - 720;
       integer_intervals_eh3 = [integer_intervals_eh3, interval_eh3];

   % 고궤도 satellite 4 위치 재생성
       interval_eh4 = interval_eh1 + random_number_3; 
       interval_eh4(interval_eh4 > 719) = interval_eh4(interval_eh4 > 719) - 720;
       integer_intervals_eh4 = [integer_intervals_eh4, interval_eh4];

   % 고궤도 satellite 5 위치 재생성
       interval_eh5 = interval_eh1 + random_number_4; 
       interval_eh5(interval_eh5 > 719) = interval_eh5(interval_eh5 > 719) - 720;
       integer_intervals_eh5 = [integer_intervals_eh5, interval_eh5];

   % 고궤도 satellite 6 위치 재생성
       interval_eh6 = interval_eh1 + random_number_5; 
       interval_eh6(interval_eh6 > 719) = interval_eh6(interval_eh6 > 719) - 720;
       integer_intervals_eh6 = [integer_intervals_eh6, interval_eh6];

   % 고궤도 satellite 7 위치 재생성
       interval_eh7 = interval_eh1 + random_number_6; 
       interval_eh7(interval_eh7 > 719) = interval_eh7(interval_eh7 > 719) - 720;
       integer_intervals_eh7 = [integer_intervals_eh7, interval_eh7];
   end


   sum_intervals_eh = ...
       [integer_intervals_eh1, ...
       integer_intervals_eh2, ...
       integer_intervals_eh3, ...
       integer_intervals_eh4, ...
       integer_intervals_eh5, ...
       integer_intervals_eh6, ...
       integer_intervals_eh7, ...
       ];


unique_intervals_eh = unique(sort(sum_intervals_eh));

Cov_eh7 = length(unique_intervals_eh) / 720 * 100;

%%
%% K_Satellite Setting(L0 / H8)

for n = 0:719
integer_intervals_kh1 = [];
integer_intervals_kh2 = [];
integer_intervals_kh3 = [];
integer_intervals_kh4 = [];
integer_intervals_kh5 = [];
integer_intervals_kh6 = [];
integer_intervals_kh7 = [];
integer_intervals_kh8 = [];

  
   for i = 1:length(ff_0)
   % 고궤도 satellite 1 위치 생성
       interval_kh1 = start_points1(i) :end_points1(i) ;
       interval_kh1(interval_kh1 > 719) = interval_kh1(interval_kh1 > 719) - 720; 
       integer_intervals_kh1 = [integer_intervals_kh1, interval_kh1];

   % 고궤도 satellite 2 위치 생성
       interval_kh2 = interval_kh1 + random_number_1; 
       interval_kh2(interval_kh2 > 719) = interval_kh2(interval_kh2 > 719) - 720;
       integer_intervals_kh2 = [integer_intervals_kh2, interval_kh2];

   % 고궤도 satellite 3 위치 생성
       interval_kh3 = interval_kh1 + random_number_2; 
       interval_kh3(interval_kh3 > 719) = interval_kh3(interval_kh3 > 719) - 720;
       integer_intervals_kh3 = [integer_intervals_kh3, interval_kh3];

   % 고궤도 satellite 4 위치 생성
       interval_kh4 = interval_kh1 + random_number_3; 
       interval_kh4(interval_kh4 > 719) = interval_kh4(interval_kh4 > 719) - 720;
       integer_intervals_kh4 = [integer_intervals_kh4, interval_kh4];

   % 고궤도 satellite 5 위치 생성
       interval_kh5 = interval_kh1 + random_number_4; 
       interval_kh5(interval_kh5 > 719) = interval_kh5(interval_kh5 > 719) - 720;
       integer_intervals_kh5 = [integer_intervals_kh5, interval_kh5];

   % 고궤도 satellite 6 위치 생성
       interval_kh6 = interval_kh1 + random_number_5; 
       interval_kh6(interval_kh6 > 719) = interval_kh6(interval_kh6 > 719) - 720;
       integer_intervals_kh6 = [integer_intervals_kh6, interval_kh6];

   % 고궤도 satellite 7 위치 생성
       interval_kh7 = interval_kh1 + random_number_6; 
       interval_kh7(interval_kh7 > 719) = interval_kh7(interval_kh7 > 719) - 720;
       integer_intervals_kh7 = [integer_intervals_kh7, interval_kh7];

   % 고궤도 satellite 8 랜덤 위치 생성
       interval_kh8 = interval_kh1 + n; 
       interval_kh8(interval_kh8 > 719) = interval_kh8(interval_kh8 > 719) - 720;
       integer_intervals_kh8 = [integer_intervals_kh8, interval_kh8];
   end


   sum_intervals_kh = ...
       [integer_intervals_kh1, ...
       integer_intervals_kh2, ...
       integer_intervals_kh3, ...
       integer_intervals_kh4, ...
       integer_intervals_kh5, ...
       integer_intervals_kh6, ...
       integer_intervals_kh7, ...
       integer_intervals_kh8, ...
       ];


unique_intervals_kh = unique(sort(sum_intervals_kh));

difference_kh = length(sum_intervals_kh) - length(unique_intervals_kh) ;

result_kh7 = [result_kh7; n, difference_kh];

end

%% E_Satellite Setting(L0 / H8)

for n = 0:719
integer_intervals_eh1 = [];
integer_intervals_eh2 = [];
integer_intervals_eh3 = [];
integer_intervals_eh4 = [];
integer_intervals_eh5 = [];
integer_intervals_eh6 = [];
integer_intervals_eh7 = [];
integer_intervals_eh8 = [];

  
   for i = 1:length(hh_0)
   % 고궤도 satellite 1 위치 생성
       interval_eh1 = start_points3(i) :end_points3(i) ;
       interval_eh1(interval_eh1 > 719) = interval_eh1(interval_eh1 > 719) - 720; 
       integer_intervals_eh1 = [integer_intervals_eh1, interval_eh1];

   % 고궤도 satellite 2 위치 생성
       interval_eh2 = interval_eh1 + random_number_1; 
       interval_eh2(interval_eh2 > 719) = interval_eh2(interval_eh2 > 719) - 720;
       integer_intervals_eh2 = [integer_intervals_eh2, interval_eh2];

   % 고궤도 satellite 3 위치 생성
       interval_eh3 = interval_eh1 + random_number_2; 
       interval_eh3(interval_eh3 > 719) = interval_eh3(interval_eh3 > 719) - 720;
       integer_intervals_eh3 = [integer_intervals_eh3, interval_eh3];

   % 고궤도 satellite 4 위치 생성
       interval_eh4 = interval_eh1 + random_number_3; 
       interval_eh4(interval_eh4 > 719) = interval_eh4(interval_eh4 > 719) - 720;
       integer_intervals_eh4 = [integer_intervals_eh4, interval_eh4];

   % 고궤도 satellite 5 위치 생성
       interval_eh5 = interval_eh1 + random_number_4; 
       interval_eh5(interval_eh5 > 719) = interval_eh5(interval_eh5 > 719) - 720;
       integer_intervals_eh5 = [integer_intervals_eh5, interval_eh5];

   % 고궤도 satellite 6 위치 생성
       interval_eh6 = interval_eh1 + random_number_5; 
       interval_eh6(interval_eh6 > 719) = interval_eh6(interval_eh6 > 719) - 720;
       integer_intervals_eh6 = [integer_intervals_eh6, interval_eh6];

   % 고궤도 satellite 7 위치 생성
       interval_eh7 = interval_eh1 + random_number_6; 
       interval_eh7(interval_eh7 > 719) = interval_eh7(interval_eh7 > 719) - 720;
       integer_intervals_eh7 = [integer_intervals_eh7, interval_eh7];

   % 고궤도 satellite 8 랜덤 위치 생성
       interval_eh8 = interval_eh1 + n; 
       interval_eh8(interval_eh8 > 719) = interval_eh8(interval_eh8 > 719) - 720;
       integer_intervals_eh8 = [integer_intervals_eh8, interval_eh8];
   end


   sum_intervals_eh = ...
       [integer_intervals_eh1, ...
       integer_intervals_eh2, ...
       integer_intervals_eh3, ...
       integer_intervals_eh4, ...
       integer_intervals_eh5, ...
       integer_intervals_eh6, ...
       integer_intervals_eh7, ...
       integer_intervals_eh8, ...
       ];


unique_intervals_eh = unique(sort(sum_intervals_eh));

difference_eh = length(sum_intervals_eh) - length(unique_intervals_eh) ;

result_eh7 = [result_eh7; n, difference_eh];

end

%% 공통

result_sum7 = [result_kh7(:,1), result_kh7(:,2) + result_eh7(:,2)];


rows_with_min_sum7 = result_sum7(result_sum7(:, 2) == min(result_sum7(:, 2)), :);

k_sum7 = rows_with_min_sum7(:,1);

random_index_7 = randi([1, length(k_sum7)]);

random_number_7 = k_sum7(random_index_7);

%% K_renewal 8
integer_intervals_kh1 = [];
integer_intervals_kh2 = [];
integer_intervals_kh3 = [];
integer_intervals_kh4 = [];
integer_intervals_kh5 = [];
integer_intervals_kh6 = [];
integer_intervals_kh7 = [];
integer_intervals_kh8 = [];

   for i = 1:length(ff_0)
   % 고궤도 satellite 1 위치 생성
       interval_kh1 = start_points1(i) :end_points1(i) ;
       interval_kh1(interval_kh1 > 719) = interval_kh1(interval_kh1 > 719) - 7320; 
       integer_intervals_kh1 = [integer_intervals_kh1, interval_kh1];

   % 고궤도 satellite 2 위치 재생성
       interval_kh2 = interval_kh1 + random_number_1; 
       interval_kh2(interval_kh2 > 719) = interval_kh2(interval_kh2 > 719) - 720;
       integer_intervals_kh2 = [integer_intervals_kh2, interval_kh2];

   % 고궤도 satellite 3 위치 재생성
       interval_kh3 = interval_kh1 + random_number_2; 
       interval_kh3(interval_kh3 > 719) = interval_kh3(interval_kh3 > 719) - 720;
       integer_intervals_kh3 = [integer_intervals_kh3, interval_kh3];

   % 고궤도 satellite 4 위치 재생성
       interval_kh4 = interval_kh1 + random_number_3; 
       interval_kh4(interval_kh4 > 719) = interval_kh4(interval_kh4 > 719) - 720;
       integer_intervals_kh4 = [integer_intervals_kh4, interval_kh4];

   % 고궤도 satellite 5 위치 재생성
       interval_kh5 = interval_kh1 + random_number_4; 
       interval_kh5(interval_kh5 > 719) = interval_kh5(interval_kh5 > 719) - 720;
       integer_intervals_kh5 = [integer_intervals_kh5, interval_kh5];

   % 고궤도 satellite 6 위치 재생성
       interval_kh6 = interval_kh1 + random_number_5; 
       interval_kh6(interval_kh6 > 719) = interval_kh6(interval_kh6 > 719) - 720;
       integer_intervals_kh6 = [integer_intervals_kh6, interval_kh6];

   % 고궤도 satellite 7 위치 재생성
       interval_kh7 = interval_kh1 + random_number_6; 
       interval_kh7(interval_kh7 > 719) = interval_kh7(interval_kh7 > 719) - 720;
       integer_intervals_kh7 = [integer_intervals_kh7, interval_kh7];

   % 고궤도 satellite 8 위치 재생성
       interval_kh8 = interval_kh1 + random_number_7; 
       interval_kh8(interval_kh8 > 719) = interval_kh8(interval_kh8 > 719) - 720;
       integer_intervals_kh8 = [integer_intervals_kh8, interval_kh8];
   end


   sum_intervals_kh = ...
       [integer_intervals_kh1, ...
       integer_intervals_kh2, ...
       integer_intervals_kh3, ...
       integer_intervals_kh4, ...
       integer_intervals_kh5, ...
       integer_intervals_kh6, ...
       integer_intervals_kh7, ...
       integer_intervals_kh8, ...
       ];


unique_intervals_kh = unique(sort(sum_intervals_kh));

Cov_kh8 = length(unique_intervals_kh) / 720 * 100;

%% E_renewal 8
integer_intervals_eh1 = [];
integer_intervals_eh2 = [];
integer_intervals_eh3 = [];
integer_intervals_eh4 = [];
integer_intervals_eh5 = [];
integer_intervals_eh6 = [];
integer_intervals_eh7 = [];
integer_intervals_eh8 = [];

   for i = 1:length(hh_0)
   % 고궤도 satellite 1 위치 생성
   interval_eh1 = start_points3(i) :end_points3(i) ;
       interval_eh1(interval_eh1 > 719) = interval_eh1(interval_eh1 > 719) - 720; 
       integer_intervals_eh1 = [integer_intervals_eh1, interval_eh1];

   % 고궤도 satellite 2 위치 재생성
       interval_eh2 = interval_eh1 + random_number_1; 
       interval_eh2(interval_eh2 > 719) = interval_eh2(interval_eh2 > 719) - 720;
       integer_intervals_eh2 = [integer_intervals_eh2, interval_eh2];

   % 고궤도 satellite 3 위치 재생성
       interval_eh3 = interval_eh1 + random_number_2; 
       interval_eh3(interval_eh3 > 719) = interval_eh3(interval_eh3 > 719) - 720;
       integer_intervals_eh3 = [integer_intervals_eh3, interval_eh3];

   % 고궤도 satellite 4 위치 재생성
       interval_eh4 = interval_eh1 + random_number_3; 
       interval_eh4(interval_eh4 > 719) = interval_eh4(interval_eh4 > 719) - 720;
       integer_intervals_eh4 = [integer_intervals_eh4, interval_eh4];

   % 고궤도 satellite 5 위치 재생성
       interval_eh5 = interval_eh1 + random_number_4; 
       interval_eh5(interval_eh5 > 719) = interval_eh5(interval_eh5 > 719) - 720;
       integer_intervals_eh5 = [integer_intervals_eh5, interval_eh5];

   % 고궤도 satellite 6 위치 재생성
       interval_eh6 = interval_eh1 + random_number_5; 
       interval_eh6(interval_eh6 > 719) = interval_eh6(interval_eh6 > 719) - 720;
       integer_intervals_eh6 = [integer_intervals_eh6, interval_eh6];

   % 고궤도 satellite 7 위치 재생성
       interval_eh7 = interval_eh1 + random_number_6; 
       interval_eh7(interval_eh7 > 719) = interval_eh7(interval_eh7 > 719) - 720;
       integer_intervals_eh7 = [integer_intervals_eh7, interval_eh7];

   % 고궤도 satellite 8 위치 재생성
       interval_eh8 = interval_eh1 + random_number_7; 
       interval_eh8(interval_eh8 > 719) = interval_eh8(interval_eh8 > 719) - 720;
       integer_intervals_eh8 = [integer_intervals_eh8, interval_eh8];
   end


   sum_intervals_eh = ...
       [integer_intervals_eh1, ...
       integer_intervals_eh2, ...
       integer_intervals_eh3, ...
       integer_intervals_eh4, ...
       integer_intervals_eh5, ...
       integer_intervals_eh6, ...
       integer_intervals_eh7, ...
       integer_intervals_eh8, ...
       ];


unique_intervals_eh = unique(sort(sum_intervals_eh));

Cov_eh8 = length(unique_intervals_eh) / 720 * 100;


%%
%% K_Satellite Setting(L0 / H9)

for n = 0:719
integer_intervals_kh1 = [];
integer_intervals_kh2 = [];
integer_intervals_kh3 = [];
integer_intervals_kh4 = [];
integer_intervals_kh5 = [];
integer_intervals_kh6 = [];
integer_intervals_kh7 = [];
integer_intervals_kh8 = [];
integer_intervals_kh9 = [];

  
   for i = 1:length(ff_0)
   % 고궤도 satellite 1 위치 생성
       interval_kh1 = start_points1(i) :end_points1(i) ;
       interval_kh1(interval_kh1 > 719) = interval_kh1(interval_kh1 > 719) - 720; 
       integer_intervals_kh1 = [integer_intervals_kh1, interval_kh1];

   % 고궤도 satellite 2 위치 생성
       interval_kh2 = interval_kh1 + random_number_1; 
       interval_kh2(interval_kh2 > 719) = interval_kh2(interval_kh2 > 719) - 720;
       integer_intervals_kh2 = [integer_intervals_kh2, interval_kh2];

   % 고궤도 satellite 3 위치 생성
       interval_kh3 = interval_kh1 + random_number_2; 
       interval_kh3(interval_kh3 > 719) = interval_kh3(interval_kh3 > 719) - 720;
       integer_intervals_kh3 = [integer_intervals_kh3, interval_kh3];

   % 고궤도 satellite 4 위치 생성
       interval_kh4 = interval_kh1 + random_number_3; 
       interval_kh4(interval_kh4 > 719) = interval_kh4(interval_kh4 > 719) - 720;
       integer_intervals_kh4 = [integer_intervals_kh4, interval_kh4];

   % 고궤도 satellite 5 위치 생성
       interval_kh5 = interval_kh1 + random_number_4; 
       interval_kh5(interval_kh5 > 719) = interval_kh5(interval_kh5 > 719) - 720;
       integer_intervals_kh5 = [integer_intervals_kh5, interval_kh5];

   % 고궤도 satellite 6 위치 생성
       interval_kh6 = interval_kh1 + random_number_5; 
       interval_kh6(interval_kh6 > 719) = interval_kh6(interval_kh6 > 719) - 720;
       integer_intervals_kh6 = [integer_intervals_kh6, interval_kh6];

   % 고궤도 satellite 7 위치 생성
       interval_kh7 = interval_kh1 + random_number_6; 
       interval_kh7(interval_kh7 > 719) = interval_kh7(interval_kh7 > 719) - 720;
       integer_intervals_kh7 = [integer_intervals_kh7, interval_kh7];

   % 고궤도 satellite 8 위치 생성
       interval_kh8 = interval_kh1 + random_number_7; 
       interval_kh8(interval_kh8 > 719) = interval_kh8(interval_kh8 > 719) - 720;
       integer_intervals_kh8 = [integer_intervals_kh8, interval_kh8];

   % 고궤도 satellite 9 랜덤 위치 생성
       interval_kh9 = interval_kh1 + n; 
       interval_kh9(interval_kh9 > 719) = interval_kh9(interval_kh9 > 719) - 720;
       integer_intervals_kh9 = [integer_intervals_kh9, interval_kh9];
   end


   sum_intervals_kh = ...
       [integer_intervals_kh1, ...
       integer_intervals_kh2, ...
       integer_intervals_kh3, ...
       integer_intervals_kh4, ...
       integer_intervals_kh5, ...
       integer_intervals_kh6, ...
       integer_intervals_kh7, ...
       integer_intervals_kh8, ...
       integer_intervals_kh9, ...
       ];


unique_intervals_kh = unique(sort(sum_intervals_kh));

difference_kh = length(sum_intervals_kh) - length(unique_intervals_kh) ;

result_kh8 = [result_kh8; n, difference_kh];

end

%% E_Satellite Setting(L0 / H9)

for n = 0:719
integer_intervals_eh1 = [];
integer_intervals_eh2 = [];
integer_intervals_eh3 = [];
integer_intervals_eh4 = [];
integer_intervals_eh5 = [];
integer_intervals_eh6 = [];
integer_intervals_eh7 = [];
integer_intervals_eh8 = [];
integer_intervals_eh9 = [];

  
   for i = 1:length(hh_0)
   % 고궤도 satellite 1 위치 생성
       interval_eh1 = start_points3(i) :end_points3(i) ;
       interval_eh1(interval_eh1 > 719) = interval_eh1(interval_eh1 > 719) - 720; 
       integer_intervals_eh1 = [integer_intervals_eh1, interval_eh1];

   % 고궤도 satellite 2 위치 생성
       interval_eh2 = interval_eh1 + random_number_1; 
       interval_eh2(interval_eh2 > 719) = interval_eh2(interval_eh2 > 719) - 720;
       integer_intervals_eh2 = [integer_intervals_eh2, interval_eh2];

   % 고궤도 satellite 3 위치 생성
       interval_eh3 = interval_eh1 + random_number_2; 
       interval_eh3(interval_eh3 > 719) = interval_eh3(interval_eh3 > 719) - 720;
       integer_intervals_eh3 = [integer_intervals_eh3, interval_eh3];

   % 고궤도 satellite 4 위치 생성
       interval_eh4 = interval_eh1 + random_number_3; 
       interval_eh4(interval_eh4 > 719) = interval_eh4(interval_eh4 > 719) - 720;
       integer_intervals_eh4 = [integer_intervals_eh4, interval_eh4];

   % 고궤도 satellite 5 위치 생성
       interval_eh5 = interval_eh1 + random_number_4; 
       interval_eh5(interval_eh5 > 719) = interval_eh5(interval_eh5 > 719) - 720;
       integer_intervals_eh5 = [integer_intervals_eh5, interval_eh5];

   % 고궤도 satellite 6 위치 생성
       interval_eh6 = interval_eh1 + random_number_5; 
       interval_eh6(interval_eh6 > 719) = interval_eh6(interval_eh6 > 719) - 720;
       integer_intervals_eh6 = [integer_intervals_eh6, interval_eh6];

   % 고궤도 satellite 7 위치 생성
       interval_eh7 = interval_eh1 + random_number_6; 
       interval_eh7(interval_eh7 > 719) = interval_eh7(interval_eh7 > 719) - 720;
       integer_intervals_eh7 = [integer_intervals_eh7, interval_eh7];

   % 고궤도 satellite 8 위치 생성
       interval_eh8 = interval_eh1 + random_number_7; 
       interval_eh8(interval_eh8 > 719) = interval_eh8(interval_eh8 > 719) - 720;
       integer_intervals_eh8 = [integer_intervals_eh8, interval_eh8];

   % 고궤도 satellite 9 랜덤 위치 생성
       interval_eh9 = interval_eh1 + n; 
       interval_eh9(interval_eh9 > 719) = interval_eh9(interval_eh9 > 719) - 720;
       integer_intervals_eh9 = [integer_intervals_eh9, interval_eh9];
   end


   sum_intervals_eh = ...
       [integer_intervals_eh1, ...
       integer_intervals_eh2, ...
       integer_intervals_eh3, ...
       integer_intervals_eh4, ...
       integer_intervals_eh5, ...
       integer_intervals_eh6, ...
       integer_intervals_eh7, ...
       integer_intervals_eh8, ...
       integer_intervals_eh9, ...
       ];


unique_intervals_eh = unique(sort(sum_intervals_eh));

difference_eh = length(sum_intervals_eh) - length(unique_intervals_eh) ;

result_eh8 = [result_eh8; n, difference_eh];

end

%% 공통

result_sum8 = [result_kh8(:,1), result_kh8(:,2) + result_eh8(:,2)];


rows_with_min_sum8 = result_sum8(result_sum8(:, 2) == min(result_sum8(:, 2)), :);

k_sum8 = rows_with_min_sum8(:,1);

random_index_8 = randi([1, length(k_sum8)]);

random_number_8 = k_sum8(random_index_8);

%% K_renewal 9
integer_intervals_kh1 = [];
integer_intervals_kh2 = [];
integer_intervals_kh3 = [];
integer_intervals_kh4 = [];
integer_intervals_kh5 = [];
integer_intervals_kh6 = [];
integer_intervals_kh7 = [];
integer_intervals_kh8 = [];
integer_intervals_kh9 = [];

   for i = 1:length(ff_0)
   % 고궤도 satellite 1 위치 생성
       interval_kh1 = start_points1(i) :end_points1(i) ;
       interval_kh1(interval_kh1 > 719) = interval_kh1(interval_kh1 > 719) - 7320; 
       integer_intervals_kh1 = [integer_intervals_kh1, interval_kh1];

   % 고궤도 satellite 2 위치 재생성
       interval_kh2 = interval_kh1 + random_number_1; 
       interval_kh2(interval_kh2 > 719) = interval_kh2(interval_kh2 > 719) - 720;
       integer_intervals_kh2 = [integer_intervals_kh2, interval_kh2];

   % 고궤도 satellite 3 위치 재생성
       interval_kh3 = interval_kh1 + random_number_2; 
       interval_kh3(interval_kh3 > 719) = interval_kh3(interval_kh3 > 719) - 720;
       integer_intervals_kh3 = [integer_intervals_kh3, interval_kh3];

   % 고궤도 satellite 4 위치 재생성
       interval_kh4 = interval_kh1 + random_number_3; 
       interval_kh4(interval_kh4 > 719) = interval_kh4(interval_kh4 > 719) - 720;
       integer_intervals_kh4 = [integer_intervals_kh4, interval_kh4];

   % 고궤도 satellite 5 위치 재생성
       interval_kh5 = interval_kh1 + random_number_4; 
       interval_kh5(interval_kh5 > 719) = interval_kh5(interval_kh5 > 719) - 720;
       integer_intervals_kh5 = [integer_intervals_kh5, interval_kh5];

   % 고궤도 satellite 6 위치 재생성
       interval_kh6 = interval_kh1 + random_number_5; 
       interval_kh6(interval_kh6 > 719) = interval_kh6(interval_kh6 > 719) - 720;
       integer_intervals_kh6 = [integer_intervals_kh6, interval_kh6];

   % 고궤도 satellite 7 위치 재생성
       interval_kh7 = interval_kh1 + random_number_6; 
       interval_kh7(interval_kh7 > 719) = interval_kh7(interval_kh7 > 719) - 720;
       integer_intervals_kh7 = [integer_intervals_kh7, interval_kh7];

   % 고궤도 satellite 8 위치 재생성
       interval_kh8 = interval_kh1 + random_number_7; 
       interval_kh8(interval_kh8 > 719) = interval_kh8(interval_kh8 > 719) - 720;
       integer_intervals_kh8 = [integer_intervals_kh8, interval_kh8];

   % 고궤도 satellite 9 위치 재생성
       interval_kh9 = interval_kh1 + random_number_8; 
       interval_kh9(interval_kh9 > 719) = interval_kh9(interval_kh9 > 719) - 720;
       integer_intervals_kh9 = [integer_intervals_kh9, interval_kh9];
   end


   sum_intervals_kh = ...
       [integer_intervals_kh1, ...
       integer_intervals_kh2, ...
       integer_intervals_kh3, ...
       integer_intervals_kh4, ...
       integer_intervals_kh5, ...
       integer_intervals_kh6, ...
       integer_intervals_kh7, ...
       integer_intervals_kh8, ...
       integer_intervals_kh9, ...
       ];


unique_intervals_kh = unique(sort(sum_intervals_kh));

Cov_kh9 = length(unique_intervals_kh) / 720 * 100;

%% E_renewal 9
integer_intervals_eh1 = [];
integer_intervals_eh2 = [];
integer_intervals_eh3 = [];
integer_intervals_eh4 = [];
integer_intervals_eh5 = [];
integer_intervals_eh6 = [];
integer_intervals_eh7 = [];
integer_intervals_eh8 = [];
integer_intervals_eh9 = [];

   for i = 1:length(hh_0)
   % 고궤도 satellite 1 위치 생성
   interval_eh1 = start_points3(i) :end_points3(i) ;
       interval_eh1(interval_eh1 > 719) = interval_eh1(interval_eh1 > 719) - 720; 
       integer_intervals_eh1 = [integer_intervals_eh1, interval_eh1];

   % 고궤도 satellite 2 위치 재생성
       interval_eh2 = interval_eh1 + random_number_1; 
       interval_eh2(interval_eh2 > 719) = interval_eh2(interval_eh2 > 719) - 720;
       integer_intervals_eh2 = [integer_intervals_eh2, interval_eh2];

   % 고궤도 satellite 3 위치 재생성
       interval_eh3 = interval_eh1 + random_number_2; 
       interval_eh3(interval_eh3 > 719) = interval_eh3(interval_eh3 > 719) - 720;
       integer_intervals_eh3 = [integer_intervals_eh3, interval_eh3];

   % 고궤도 satellite 4 위치 재생성
       interval_eh4 = interval_eh1 + random_number_3; 
       interval_eh4(interval_eh4 > 719) = interval_eh4(interval_eh4 > 719) - 720;
       integer_intervals_eh4 = [integer_intervals_eh4, interval_eh4];

   % 고궤도 satellite 5 위치 재생성
       interval_eh5 = interval_eh1 + random_number_4; 
       interval_eh5(interval_eh5 > 719) = interval_eh5(interval_eh5 > 719) - 720;
       integer_intervals_eh5 = [integer_intervals_eh5, interval_eh5];

   % 고궤도 satellite 6 위치 재생성
       interval_eh6 = interval_eh1 + random_number_5; 
       interval_eh6(interval_eh6 > 719) = interval_eh6(interval_eh6 > 719) - 720;
       integer_intervals_eh6 = [integer_intervals_eh6, interval_eh6];

   % 고궤도 satellite 7 위치 재생성
       interval_eh7 = interval_eh1 + random_number_6; 
       interval_eh7(interval_eh7 > 719) = interval_eh7(interval_eh7 > 719) - 720;
       integer_intervals_eh7 = [integer_intervals_eh7, interval_eh7];

   % 고궤도 satellite 8 위치 재생성
       interval_eh8 = interval_eh1 + random_number_7; 
       interval_eh8(interval_eh8 > 719) = interval_eh8(interval_eh8 > 719) - 720;
       integer_intervals_eh8 = [integer_intervals_eh8, interval_eh8];

   % 고궤도 satellite 9 위치 재생성
       interval_eh9 = interval_eh1 + random_number_8; 
       interval_eh9(interval_eh9 > 719) = interval_eh9(interval_eh9 > 719) - 720;
       integer_intervals_eh9 = [integer_intervals_eh9, interval_eh9];
   end
   
   sum_intervals_eh = ...
       [integer_intervals_eh1, ...
       integer_intervals_eh2, ...
       integer_intervals_eh3, ...
       integer_intervals_eh4, ...
       integer_intervals_eh5, ...
       integer_intervals_eh6, ...
       integer_intervals_eh7, ...
       integer_intervals_eh8, ...
       integer_intervals_eh9, ...
       ];


unique_intervals_eh = unique(sort(sum_intervals_eh));

Cov_eh9 = length(unique_intervals_eh) / 720 * 100;

%%
%% K_Satellite Setting(L0 / H10)

for n = 0:719
integer_intervals_kh1 = [];
integer_intervals_kh2 = [];
integer_intervals_kh3 = [];
integer_intervals_kh4 = [];
integer_intervals_kh5 = [];
integer_intervals_kh6 = [];
integer_intervals_kh7 = [];
integer_intervals_kh8 = [];
integer_intervals_kh9 = [];
integer_intervals_kh10 = [];

  
   for i = 1:length(ff_0)
   % 고궤도 satellite 1 위치 생성
       interval_kh1 = start_points1(i) :end_points1(i) ;
       interval_kh1(interval_kh1 > 719) = interval_kh1(interval_kh1 > 719) - 720; 
       integer_intervals_kh1 = [integer_intervals_kh1, interval_kh1];

   % 고궤도 satellite 2 위치 생성
       interval_kh2 = interval_kh1 + random_number_1; 
       interval_kh2(interval_kh2 > 719) = interval_kh2(interval_kh2 > 719) - 720;
       integer_intervals_kh2 = [integer_intervals_kh2, interval_kh2];

   % 고궤도 satellite 3 위치 생성
       interval_kh3 = interval_kh1 + random_number_2; 
       interval_kh3(interval_kh3 > 719) = interval_kh3(interval_kh3 > 719) - 720;
       integer_intervals_kh3 = [integer_intervals_kh3, interval_kh3];

   % 고궤도 satellite 4 위치 생성
       interval_kh4 = interval_kh1 + random_number_3; 
       interval_kh4(interval_kh4 > 719) = interval_kh4(interval_kh4 > 719) - 720;
       integer_intervals_kh4 = [integer_intervals_kh4, interval_kh4];

   % 고궤도 satellite 5 위치 생성
       interval_kh5 = interval_kh1 + random_number_4; 
       interval_kh5(interval_kh5 > 719) = interval_kh5(interval_kh5 > 719) - 720;
       integer_intervals_kh5 = [integer_intervals_kh5, interval_kh5];

   % 고궤도 satellite 6 위치 생성
       interval_kh6 = interval_kh1 + random_number_5; 
       interval_kh6(interval_kh6 > 719) = interval_kh6(interval_kh6 > 719) - 720;
       integer_intervals_kh6 = [integer_intervals_kh6, interval_kh6];

   % 고궤도 satellite 7 위치 생성
       interval_kh7 = interval_kh1 + random_number_6; 
       interval_kh7(interval_kh7 > 719) = interval_kh7(interval_kh7 > 719) - 720;
       integer_intervals_kh7 = [integer_intervals_kh7, interval_kh7];

   % 고궤도 satellite 8 위치 생성
       interval_kh8 = interval_kh1 + random_number_7; 
       interval_kh8(interval_kh8 > 719) = interval_kh8(interval_kh8 > 719) - 720;
       integer_intervals_kh8 = [integer_intervals_kh8, interval_kh8];

   % 고궤도 satellite 9 위치 생성
       interval_kh9 = interval_kh1 + random_number_8; 
       interval_kh9(interval_kh9 > 719) = interval_kh9(interval_kh9 > 719) - 720;
       integer_intervals_kh9 = [integer_intervals_kh9, interval_kh9];

   % 고궤도 satellite 10 랜덤 위치 생성
       interval_kh10 = interval_kh1 + n; 
       interval_kh10(interval_kh10 > 719) = interval_kh10(interval_kh10 > 719) - 720;
       integer_intervals_kh10 = [integer_intervals_kh10, interval_kh10];
   end


   sum_intervals_kh = ...
       [integer_intervals_kh1, ...
       integer_intervals_kh2, ...
       integer_intervals_kh3, ...
       integer_intervals_kh4, ...
       integer_intervals_kh5, ...
       integer_intervals_kh6, ...
       integer_intervals_kh7, ...
       integer_intervals_kh8, ...
       integer_intervals_kh9, ...
       integer_intervals_kh10, ...
       ];


unique_intervals_kh = unique(sort(sum_intervals_kh));

difference_kh = length(sum_intervals_kh) - length(unique_intervals_kh) ;

result_kh9 = [result_kh9; n, difference_kh];

end

%% E_Satellite Setting(L0 / H10)

for n = 0:719
integer_intervals_eh1 = [];
integer_intervals_eh2 = [];
integer_intervals_eh3 = [];
integer_intervals_eh4 = [];
integer_intervals_eh5 = [];
integer_intervals_eh6 = [];
integer_intervals_eh7 = [];
integer_intervals_eh8 = [];
integer_intervals_eh9 = [];
integer_intervals_eh10 = [];

  
   for i = 1:length(hh_0)
   % 고궤도 satellite 1 위치 생성
       interval_eh1 = start_points3(i) :end_points3(i) ;
       interval_eh1(interval_eh1 > 719) = interval_eh1(interval_eh1 > 719) - 720; 
       integer_intervals_eh1 = [integer_intervals_eh1, interval_eh1];

   % 고궤도 satellite 2 위치 생성
       interval_eh2 = interval_eh1 + random_number_1; 
       interval_eh2(interval_eh2 > 719) = interval_eh2(interval_eh2 > 719) - 720;
       integer_intervals_eh2 = [integer_intervals_eh2, interval_eh2];

   % 고궤도 satellite 3 위치 생성
       interval_eh3 = interval_eh1 + random_number_2; 
       interval_eh3(interval_eh3 > 719) = interval_eh3(interval_eh3 > 719) - 720;
       integer_intervals_eh3 = [integer_intervals_eh3, interval_eh3];

   % 고궤도 satellite 4 위치 생성
       interval_eh4 = interval_eh1 + random_number_3; 
       interval_eh4(interval_eh4 > 719) = interval_eh4(interval_eh4 > 719) - 720;
       integer_intervals_eh4 = [integer_intervals_eh4, interval_eh4];

   % 고궤도 satellite 5 위치 생성
       interval_eh5 = interval_eh1 + random_number_4; 
       interval_eh5(interval_eh5 > 719) = interval_eh5(interval_eh5 > 719) - 720;
       integer_intervals_eh5 = [integer_intervals_eh5, interval_eh5];

   % 고궤도 satellite 6 위치 생성
       interval_eh6 = interval_eh1 + random_number_5; 
       interval_eh6(interval_eh6 > 719) = interval_eh6(interval_eh6 > 719) - 720;
       integer_intervals_eh6 = [integer_intervals_eh6, interval_eh6];

   % 고궤도 satellite 7 위치 생성
       interval_eh7 = interval_eh1 + random_number_6; 
       interval_eh7(interval_eh7 > 719) = interval_eh7(interval_eh7 > 719) - 720;
       integer_intervals_eh7 = [integer_intervals_eh7, interval_eh7];

   % 고궤도 satellite 8 위치 생성
       interval_eh8 = interval_eh1 + random_number_7; 
       interval_eh8(interval_eh8 > 719) = interval_eh8(interval_eh8 > 719) - 720;
       integer_intervals_eh8 = [integer_intervals_eh8, interval_eh8];

   % 고궤도 satellite 9 위치 생성
       interval_eh9 = interval_eh1 + random_number_8; 
       interval_eh9(interval_eh9 > 719) = interval_eh9(interval_eh9 > 719) - 720;
       integer_intervals_eh9 = [integer_intervals_eh9, interval_eh9];

   % 고궤도 satellite 10 랜덤 위치 생성
       interval_eh10 = interval_eh1 + n; 
       interval_eh10(interval_eh10 > 719) = interval_eh10(interval_eh10 > 719) - 720;
       integer_intervals_eh10 = [integer_intervals_eh10, interval_eh10];
   end


   sum_intervals_eh = ...
       [integer_intervals_eh1, ...
       integer_intervals_eh2, ...
       integer_intervals_eh3, ...
       integer_intervals_eh4, ...
       integer_intervals_eh5, ...
       integer_intervals_eh6, ...
       integer_intervals_eh7, ...
       integer_intervals_eh8, ...
       integer_intervals_eh9, ...
       integer_intervals_eh10, ...
       ];


unique_intervals_eh = unique(sort(sum_intervals_eh));

difference_eh = length(sum_intervals_eh) - length(unique_intervals_eh) ;

result_eh9 = [result_eh9; n, difference_eh];

end

%% 공통

result_sum9 = [result_kh9(:,1), result_kh9(:,2) + result_eh9(:,2)];


rows_with_min_sum9 = result_sum9(result_sum9(:, 2) == min(result_sum9(:, 2)), :);

k_sum9 = rows_with_min_sum9(:,1);

random_index_9 = randi([1, length(k_sum9)]);

random_number_9 = k_sum9(random_index_9);

%% K_renewal 10
integer_intervals_kh1 = [];
integer_intervals_kh2 = [];
integer_intervals_kh3 = [];
integer_intervals_kh4 = [];
integer_intervals_kh5 = [];
integer_intervals_kh6 = [];
integer_intervals_kh7 = [];
integer_intervals_kh8 = [];
integer_intervals_kh9 = [];
integer_intervals_kh10 = [];

   for i = 1:length(ff_0)
   % 고궤도 satellite 1 위치 생성
       interval_kh1 = start_points1(i) :end_points1(i) ;
       interval_kh1(interval_kh1 > 719) = interval_kh1(interval_kh1 > 719) - 7320; 
       integer_intervals_kh1 = [integer_intervals_kh1, interval_kh1];

   % 고궤도 satellite 2 위치 재생성
       interval_kh2 = interval_kh1 + random_number_1; 
       interval_kh2(interval_kh2 > 719) = interval_kh2(interval_kh2 > 719) - 720;
       integer_intervals_kh2 = [integer_intervals_kh2, interval_kh2];

   % 고궤도 satellite 3 위치 재생성
       interval_kh3 = interval_kh1 + random_number_2; 
       interval_kh3(interval_kh3 > 719) = interval_kh3(interval_kh3 > 719) - 720;
       integer_intervals_kh3 = [integer_intervals_kh3, interval_kh3];

   % 고궤도 satellite 4 위치 재생성
       interval_kh4 = interval_kh1 + random_number_3; 
       interval_kh4(interval_kh4 > 719) = interval_kh4(interval_kh4 > 719) - 720;
       integer_intervals_kh4 = [integer_intervals_kh4, interval_kh4];

   % 고궤도 satellite 5 위치 재생성
       interval_kh5 = interval_kh1 + random_number_4; 
       interval_kh5(interval_kh5 > 719) = interval_kh5(interval_kh5 > 719) - 720;
       integer_intervals_kh5 = [integer_intervals_kh5, interval_kh5];

   % 고궤도 satellite 6 위치 재생성
       interval_kh6 = interval_kh1 + random_number_5; 
       interval_kh6(interval_kh6 > 719) = interval_kh6(interval_kh6 > 719) - 720;
       integer_intervals_kh6 = [integer_intervals_kh6, interval_kh6];

   % 고궤도 satellite 7 위치 재생성
       interval_kh7 = interval_kh1 + random_number_6; 
       interval_kh7(interval_kh7 > 719) = interval_kh7(interval_kh7 > 719) - 720;
       integer_intervals_kh7 = [integer_intervals_kh7, interval_kh7];

   % 고궤도 satellite 8 위치 재생성
       interval_kh8 = interval_kh1 + random_number_7; 
       interval_kh8(interval_kh8 > 719) = interval_kh8(interval_kh8 > 719) - 720;
       integer_intervals_kh8 = [integer_intervals_kh8, interval_kh8];

   % 고궤도 satellite 9 위치 재생성
       interval_kh9 = interval_kh1 + random_number_8; 
       interval_kh9(interval_kh9 > 719) = interval_kh9(interval_kh9 > 719) - 720;
       integer_intervals_kh9 = [integer_intervals_kh9, interval_kh9];

   % 고궤도 satellite 10 위치 재생성
       interval_kh10 = interval_kh1 + random_number_9; 
       interval_kh10(interval_kh10 > 719) = interval_kh10(interval_kh10 > 719) - 720;
       integer_intervals_kh10 = [integer_intervals_kh10, interval_kh10];
   end


   sum_intervals_kh = ...
       [integer_intervals_kh1, ...
       integer_intervals_kh2, ...
       integer_intervals_kh3, ...
       integer_intervals_kh4, ...
       integer_intervals_kh5, ...
       integer_intervals_kh6, ...
       integer_intervals_kh7, ...
       integer_intervals_kh8, ...
       integer_intervals_kh9, ...
       integer_intervals_kh10, ...
       ];


unique_intervals_kh = unique(sort(sum_intervals_kh));

Cov_kh10 = length(unique_intervals_kh) / 720 * 100;

%% E_renewal 10
integer_intervals_eh1 = [];
integer_intervals_eh2 = [];
integer_intervals_eh3 = [];
integer_intervals_eh4 = [];
integer_intervals_eh5 = [];
integer_intervals_eh6 = [];
integer_intervals_eh7 = [];
integer_intervals_eh8 = [];
integer_intervals_eh9 = [];
integer_intervals_eh10 = [];

   for i = 1:length(hh_0)
   % 고궤도 satellite 1 위치 생성
   interval_eh1 = start_points3(i) :end_points3(i) ;
       interval_eh1(interval_eh1 > 719) = interval_eh1(interval_eh1 > 719) - 720; 
       integer_intervals_eh1 = [integer_intervals_eh1, interval_eh1];

   % 고궤도 satellite 2 위치 재생성
       interval_eh2 = interval_eh1 + random_number_1; 
       interval_eh2(interval_eh2 > 719) = interval_eh2(interval_eh2 > 719) - 720;
       integer_intervals_eh2 = [integer_intervals_eh2, interval_eh2];

   % 고궤도 satellite 3 위치 재생성
       interval_eh3 = interval_eh1 + random_number_2; 
       interval_eh3(interval_eh3 > 719) = interval_eh3(interval_eh3 > 719) - 720;
       integer_intervals_eh3 = [integer_intervals_eh3, interval_eh3];

   % 고궤도 satellite 4 위치 재생성
       interval_eh4 = interval_eh1 + random_number_3; 
       interval_eh4(interval_eh4 > 719) = interval_eh4(interval_eh4 > 719) - 720;
       integer_intervals_eh4 = [integer_intervals_eh4, interval_eh4];

   % 고궤도 satellite 5 위치 재생성
       interval_eh5 = interval_eh1 + random_number_4; 
       interval_eh5(interval_eh5 > 719) = interval_eh5(interval_eh5 > 719) - 720;
       integer_intervals_eh5 = [integer_intervals_eh5, interval_eh5];

   % 고궤도 satellite 6 위치 재생성
       interval_eh6 = interval_eh1 + random_number_5; 
       interval_eh6(interval_eh6 > 719) = interval_eh6(interval_eh6 > 719) - 720;
       integer_intervals_eh6 = [integer_intervals_eh6, interval_eh6];

   % 고궤도 satellite 7 위치 재생성
       interval_eh7 = interval_eh1 + random_number_6; 
       interval_eh7(interval_eh7 > 719) = interval_eh7(interval_eh7 > 719) - 720;
       integer_intervals_eh7 = [integer_intervals_eh7, interval_eh7];

   % 고궤도 satellite 8 위치 재생성
       interval_eh8 = interval_eh1 + random_number_7; 
       interval_eh8(interval_eh8 > 719) = interval_eh8(interval_eh8 > 719) - 720;
       integer_intervals_eh8 = [integer_intervals_eh8, interval_eh8];

   % 고궤도 satellite 9 위치 재생성
       interval_eh9 = interval_eh1 + random_number_8; 
       interval_eh9(interval_eh9 > 719) = interval_eh9(interval_eh9 > 719) - 720;
       integer_intervals_eh9 = [integer_intervals_eh9, interval_eh9];

   % 고궤도 satellite 10 위치 재생성
       interval_eh10 = interval_eh1 + random_number_9; 
       interval_eh10(interval_eh10 > 719) = interval_eh10(interval_eh10 > 719) - 720;
       integer_intervals_eh10 = [integer_intervals_eh10, interval_eh10];
   end
   
   sum_intervals_eh = ...
       [integer_intervals_eh1, ...
       integer_intervals_eh2, ...
       integer_intervals_eh3, ...
       integer_intervals_eh4, ...
       integer_intervals_eh5, ...
       integer_intervals_eh6, ...
       integer_intervals_eh7, ...
       integer_intervals_eh8, ...
       integer_intervals_eh9, ...
       integer_intervals_eh10, ...
       ];


unique_intervals_eh = unique(sort(sum_intervals_eh));

Cov_eh10 = length(unique_intervals_eh) / 720 * 100;


%%
%% K_Satellite Setting(L0 / H11)

for n = 0:719
integer_intervals_kh1 = [];
integer_intervals_kh2 = [];
integer_intervals_kh3 = [];
integer_intervals_kh4 = [];
integer_intervals_kh5 = [];
integer_intervals_kh6 = [];
integer_intervals_kh7 = [];
integer_intervals_kh8 = [];
integer_intervals_kh9 = [];
integer_intervals_kh10 = [];
integer_intervals_kh11 = [];

  
   for i = 1:length(ff_0)
   % 고궤도 satellite 1 위치 생성
       interval_kh1 = start_points1(i) :end_points1(i) ;
       interval_kh1(interval_kh1 > 719) = interval_kh1(interval_kh1 > 719) - 720; 
       integer_intervals_kh1 = [integer_intervals_kh1, interval_kh1];

   % 고궤도 satellite 2 위치 생성
       interval_kh2 = interval_kh1 + random_number_1; 
       interval_kh2(interval_kh2 > 719) = interval_kh2(interval_kh2 > 719) - 720;
       integer_intervals_kh2 = [integer_intervals_kh2, interval_kh2];

   % 고궤도 satellite 3 위치 생성
       interval_kh3 = interval_kh1 + random_number_2; 
       interval_kh3(interval_kh3 > 719) = interval_kh3(interval_kh3 > 719) - 720;
       integer_intervals_kh3 = [integer_intervals_kh3, interval_kh3];

   % 고궤도 satellite 4 위치 생성
       interval_kh4 = interval_kh1 + random_number_3; 
       interval_kh4(interval_kh4 > 719) = interval_kh4(interval_kh4 > 719) - 720;
       integer_intervals_kh4 = [integer_intervals_kh4, interval_kh4];

   % 고궤도 satellite 5 위치 생성
       interval_kh5 = interval_kh1 + random_number_4; 
       interval_kh5(interval_kh5 > 719) = interval_kh5(interval_kh5 > 719) - 720;
       integer_intervals_kh5 = [integer_intervals_kh5, interval_kh5];

   % 고궤도 satellite 6 위치 생성
       interval_kh6 = interval_kh1 + random_number_5; 
       interval_kh6(interval_kh6 > 719) = interval_kh6(interval_kh6 > 719) - 720;
       integer_intervals_kh6 = [integer_intervals_kh6, interval_kh6];

   % 고궤도 satellite 7 위치 생성
       interval_kh7 = interval_kh1 + random_number_6; 
       interval_kh7(interval_kh7 > 719) = interval_kh7(interval_kh7 > 719) - 720;
       integer_intervals_kh7 = [integer_intervals_kh7, interval_kh7];

   % 고궤도 satellite 8 위치 생성
       interval_kh8 = interval_kh1 + random_number_7; 
       interval_kh8(interval_kh8 > 719) = interval_kh8(interval_kh8 > 719) - 720;
       integer_intervals_kh8 = [integer_intervals_kh8, interval_kh8];

   % 고궤도 satellite 9 위치 생성
       interval_kh9 = interval_kh1 + random_number_8; 
       interval_kh9(interval_kh9 > 719) = interval_kh9(interval_kh9 > 719) - 720;
       integer_intervals_kh9 = [integer_intervals_kh9, interval_kh9];

   % 고궤도 satellite 10 위치 생성
       interval_kh10 = interval_kh1 + random_number_9; 
       interval_kh10(interval_kh10 > 719) = interval_kh10(interval_kh10 > 719) - 720;
       integer_intervals_kh10 = [integer_intervals_kh10, interval_kh10];

   % 고궤도 satellite 11 랜덤 위치 생성
       interval_kh11 = interval_kh1 + n; 
       interval_kh11(interval_kh11 > 719) = interval_kh11(interval_kh11 > 719) - 720;
       integer_intervals_kh11 = [integer_intervals_kh11, interval_kh11];
   end


   sum_intervals_kh = ...
       [integer_intervals_kh1, ...
       integer_intervals_kh2, ...
       integer_intervals_kh3, ...
       integer_intervals_kh4, ...
       integer_intervals_kh5, ...
       integer_intervals_kh6, ...
       integer_intervals_kh7, ...
       integer_intervals_kh8, ...
       integer_intervals_kh9, ...
       integer_intervals_kh10, ...
       integer_intervals_kh11, ...
       ];


unique_intervals_kh = unique(sort(sum_intervals_kh));

difference_kh = length(sum_intervals_kh) - length(unique_intervals_kh) ;

result_kh10 = [result_kh10; n, difference_kh];

end

%% E_Satellite Setting(L0 / H11)

for n = 0:719
integer_intervals_eh1 = [];
integer_intervals_eh2 = [];
integer_intervals_eh3 = [];
integer_intervals_eh4 = [];
integer_intervals_eh5 = [];
integer_intervals_eh6 = [];
integer_intervals_eh7 = [];
integer_intervals_eh8 = [];
integer_intervals_eh9 = [];
integer_intervals_eh10 = [];
integer_intervals_eh11 = [];

  
   for i = 1:length(hh_0)
   % 고궤도 satellite 1 위치 생성
       interval_eh1 = start_points3(i) :end_points3(i) ;
       interval_eh1(interval_eh1 > 719) = interval_eh1(interval_eh1 > 719) - 720; 
       integer_intervals_eh1 = [integer_intervals_eh1, interval_eh1];

   % 고궤도 satellite 2 위치 생성
       interval_eh2 = interval_eh1 + random_number_1; 
       interval_eh2(interval_eh2 > 719) = interval_eh2(interval_eh2 > 719) - 720;
       integer_intervals_eh2 = [integer_intervals_eh2, interval_eh2];

   % 고궤도 satellite 3 위치 생성
       interval_eh3 = interval_eh1 + random_number_2; 
       interval_eh3(interval_eh3 > 719) = interval_eh3(interval_eh3 > 719) - 720;
       integer_intervals_eh3 = [integer_intervals_eh3, interval_eh3];

   % 고궤도 satellite 4 위치 생성
       interval_eh4 = interval_eh1 + random_number_3; 
       interval_eh4(interval_eh4 > 719) = interval_eh4(interval_eh4 > 719) - 720;
       integer_intervals_eh4 = [integer_intervals_eh4, interval_eh4];

   % 고궤도 satellite 5 위치 생성
       interval_eh5 = interval_eh1 + random_number_4; 
       interval_eh5(interval_eh5 > 719) = interval_eh5(interval_eh5 > 719) - 720;
       integer_intervals_eh5 = [integer_intervals_eh5, interval_eh5];

   % 고궤도 satellite 6 위치 생성
       interval_eh6 = interval_eh1 + random_number_5; 
       interval_eh6(interval_eh6 > 719) = interval_eh6(interval_eh6 > 719) - 720;
       integer_intervals_eh6 = [integer_intervals_eh6, interval_eh6];

   % 고궤도 satellite 7 위치 생성
       interval_eh7 = interval_eh1 + random_number_6; 
       interval_eh7(interval_eh7 > 719) = interval_eh7(interval_eh7 > 719) - 720;
       integer_intervals_eh7 = [integer_intervals_eh7, interval_eh7];

   % 고궤도 satellite 8 위치 생성
       interval_eh8 = interval_eh1 + random_number_7; 
       interval_eh8(interval_eh8 > 719) = interval_eh8(interval_eh8 > 719) - 720;
       integer_intervals_eh8 = [integer_intervals_eh8, interval_eh8];

   % 고궤도 satellite 9 위치 생성
       interval_eh9 = interval_eh1 + random_number_8; 
       interval_eh9(interval_eh9 > 719) = interval_eh9(interval_eh9 > 719) - 720;
       integer_intervals_eh9 = [integer_intervals_eh9, interval_eh9];

   % 고궤도 satellite 10 위치 생성
       interval_eh10 = interval_eh1 + random_number_9; 
       interval_eh10(interval_eh10 > 719) = interval_eh10(interval_eh10 > 719) - 720;
       integer_intervals_eh10 = [integer_intervals_eh10, interval_eh10];

   % 고궤도 satellite 11 랜덤 위치 생성
       interval_eh11 = interval_eh1 + n; 
       interval_eh11(interval_eh11 > 719) = interval_eh11(interval_eh11 > 719) - 720;
       integer_intervals_eh11 = [integer_intervals_eh11, interval_eh11];
   end


   sum_intervals_eh = ...
       [integer_intervals_eh1, ...
       integer_intervals_eh2, ...
       integer_intervals_eh3, ...
       integer_intervals_eh4, ...
       integer_intervals_eh5, ...
       integer_intervals_eh6, ...
       integer_intervals_eh7, ...
       integer_intervals_eh8, ...
       integer_intervals_eh9, ...
       integer_intervals_eh10, ...
       integer_intervals_eh11, ...
       ];


unique_intervals_eh = unique(sort(sum_intervals_eh));

difference_eh = length(sum_intervals_eh) - length(unique_intervals_eh) ;

result_eh10 = [result_eh10; n, difference_eh];

end

%% 공통

result_sum10 = [result_kh10(:,1), result_kh10(:,2) + result_eh10(:,2)];


rows_with_min_sum10 = result_sum10(result_sum10(:, 2) == min(result_sum10(:, 2)), :);

k_sum10 = rows_with_min_sum10(:,1);

random_index_10 = randi([1, length(k_sum10)]);

random_number_10 = k_sum10(random_index_10);

%% K_renewal 11
integer_intervals_kh1 = [];
integer_intervals_kh2 = [];
integer_intervals_kh3 = [];
integer_intervals_kh4 = [];
integer_intervals_kh5 = [];
integer_intervals_kh6 = [];
integer_intervals_kh7 = [];
integer_intervals_kh8 = [];
integer_intervals_kh9 = [];
integer_intervals_kh10 = [];
integer_intervals_kh11 = [];

   for i = 1:length(ff_0)
   % 고궤도 satellite 1 위치 생성
       interval_kh1 = start_points1(i) :end_points1(i) ;
       interval_kh1(interval_kh1 > 719) = interval_kh1(interval_kh1 > 719) - 7320; 
       integer_intervals_kh1 = [integer_intervals_kh1, interval_kh1];

   % 고궤도 satellite 2 위치 재생성
       interval_kh2 = interval_kh1 + random_number_1; 
       interval_kh2(interval_kh2 > 719) = interval_kh2(interval_kh2 > 719) - 720;
       integer_intervals_kh2 = [integer_intervals_kh2, interval_kh2];

   % 고궤도 satellite 3 위치 재생성
       interval_kh3 = interval_kh1 + random_number_2; 
       interval_kh3(interval_kh3 > 719) = interval_kh3(interval_kh3 > 719) - 720;
       integer_intervals_kh3 = [integer_intervals_kh3, interval_kh3];

   % 고궤도 satellite 4 위치 재생성
       interval_kh4 = interval_kh1 + random_number_3; 
       interval_kh4(interval_kh4 > 719) = interval_kh4(interval_kh4 > 719) - 720;
       integer_intervals_kh4 = [integer_intervals_kh4, interval_kh4];

   % 고궤도 satellite 5 위치 재생성
       interval_kh5 = interval_kh1 + random_number_4; 
       interval_kh5(interval_kh5 > 719) = interval_kh5(interval_kh5 > 719) - 720;
       integer_intervals_kh5 = [integer_intervals_kh5, interval_kh5];

   % 고궤도 satellite 6 위치 재생성
       interval_kh6 = interval_kh1 + random_number_5; 
       interval_kh6(interval_kh6 > 719) = interval_kh6(interval_kh6 > 719) - 720;
       integer_intervals_kh6 = [integer_intervals_kh6, interval_kh6];

   % 고궤도 satellite 7 위치 재생성
       interval_kh7 = interval_kh1 + random_number_6; 
       interval_kh7(interval_kh7 > 719) = interval_kh7(interval_kh7 > 719) - 720;
       integer_intervals_kh7 = [integer_intervals_kh7, interval_kh7];

   % 고궤도 satellite 8 위치 재생성
       interval_kh8 = interval_kh1 + random_number_7; 
       interval_kh8(interval_kh8 > 719) = interval_kh8(interval_kh8 > 719) - 720;
       integer_intervals_kh8 = [integer_intervals_kh8, interval_kh8];

   % 고궤도 satellite 9 위치 재생성
       interval_kh9 = interval_kh1 + random_number_8; 
       interval_kh9(interval_kh9 > 719) = interval_kh9(interval_kh9 > 719) - 720;
       integer_intervals_kh9 = [integer_intervals_kh9, interval_kh9];

   % 고궤도 satellite 10 위치 재생성
       interval_kh10 = interval_kh1 + random_number_9; 
       interval_kh10(interval_kh10 > 719) = interval_kh10(interval_kh10 > 719) - 720;
       integer_intervals_kh10 = [integer_intervals_kh10, interval_kh10];

   % 고궤도 satellite 11 위치 재생성
       interval_kh11 = interval_kh1 + random_number_10; 
       interval_kh11(interval_kh11 > 719) = interval_kh11(interval_kh11 > 719) - 720;
       integer_intervals_kh11 = [integer_intervals_kh11, interval_kh11];
   end


   sum_intervals_kh = ...
       [integer_intervals_kh1, ...
       integer_intervals_kh2, ...
       integer_intervals_kh3, ...
       integer_intervals_kh4, ...
       integer_intervals_kh5, ...
       integer_intervals_kh6, ...
       integer_intervals_kh7, ...
       integer_intervals_kh8, ...
       integer_intervals_kh9, ...
       integer_intervals_kh10, ...
       integer_intervals_kh11, ...
       ];


unique_intervals_kh = unique(sort(sum_intervals_kh));

Cov_kh11 = length(unique_intervals_kh) / 720 * 100;

%% E_renewal 11
integer_intervals_eh1 = [];
integer_intervals_eh2 = [];
integer_intervals_eh3 = [];
integer_intervals_eh4 = [];
integer_intervals_eh5 = [];
integer_intervals_eh6 = [];
integer_intervals_eh7 = [];
integer_intervals_eh8 = [];
integer_intervals_eh9 = [];
integer_intervals_eh10 = [];
integer_intervals_eh11 = [];

   for i = 1:length(hh_0)
   % 고궤도 satellite 1 위치 생성
   interval_eh1 = start_points3(i) :end_points3(i) ;
       interval_eh1(interval_eh1 > 719) = interval_eh1(interval_eh1 > 719) - 720; 
       integer_intervals_eh1 = [integer_intervals_eh1, interval_eh1];

   % 고궤도 satellite 2 위치 재생성
       interval_eh2 = interval_eh1 + random_number_1; 
       interval_eh2(interval_eh2 > 719) = interval_eh2(interval_eh2 > 719) - 720;
       integer_intervals_eh2 = [integer_intervals_eh2, interval_eh2];

   % 고궤도 satellite 3 위치 재생성
       interval_eh3 = interval_eh1 + random_number_2; 
       interval_eh3(interval_eh3 > 719) = interval_eh3(interval_eh3 > 719) - 720;
       integer_intervals_eh3 = [integer_intervals_eh3, interval_eh3];

   % 고궤도 satellite 4 위치 재생성
       interval_eh4 = interval_eh1 + random_number_3; 
       interval_eh4(interval_eh4 > 719) = interval_eh4(interval_eh4 > 719) - 720;
       integer_intervals_eh4 = [integer_intervals_eh4, interval_eh4];

   % 고궤도 satellite 5 위치 재생성
       interval_eh5 = interval_eh1 + random_number_4; 
       interval_eh5(interval_eh5 > 719) = interval_eh5(interval_eh5 > 719) - 720;
       integer_intervals_eh5 = [integer_intervals_eh5, interval_eh5];

   % 고궤도 satellite 6 위치 재생성
       interval_eh6 = interval_eh1 + random_number_5; 
       interval_eh6(interval_eh6 > 719) = interval_eh6(interval_eh6 > 719) - 720;
       integer_intervals_eh6 = [integer_intervals_eh6, interval_eh6];

   % 고궤도 satellite 7 위치 재생성
       interval_eh7 = interval_eh1 + random_number_6; 
       interval_eh7(interval_eh7 > 719) = interval_eh7(interval_eh7 > 719) - 720;
       integer_intervals_eh7 = [integer_intervals_eh7, interval_eh7];

   % 고궤도 satellite 8 위치 재생성
       interval_eh8 = interval_eh1 + random_number_7; 
       interval_eh8(interval_eh8 > 719) = interval_eh8(interval_eh8 > 719) - 720;
       integer_intervals_eh8 = [integer_intervals_eh8, interval_eh8];

   % 고궤도 satellite 9 위치 재생성
       interval_eh9 = interval_eh1 + random_number_8; 
       interval_eh9(interval_eh9 > 719) = interval_eh9(interval_eh9 > 719) - 720;
       integer_intervals_eh9 = [integer_intervals_eh9, interval_eh9];

   % 고궤도 satellite 10 위치 재생성
       interval_eh10 = interval_eh1 + random_number_9; 
       interval_eh10(interval_eh10 > 719) = interval_eh10(interval_eh10 > 719) - 720;
       integer_intervals_eh10 = [integer_intervals_eh10, interval_eh10];

   % 고궤도 satellite 11 위치 재생성
       interval_eh11 = interval_eh1 + random_number_10; 
       interval_eh11(interval_eh11 > 719) = interval_eh11(interval_eh11 > 719) - 720;
       integer_intervals_eh11 = [integer_intervals_eh11, interval_eh11];
   end
   
   sum_intervals_eh = ...
       [integer_intervals_eh1, ...
       integer_intervals_eh2, ...
       integer_intervals_eh3, ...
       integer_intervals_eh4, ...
       integer_intervals_eh5, ...
       integer_intervals_eh6, ...
       integer_intervals_eh7, ...
       integer_intervals_eh8, ...
       integer_intervals_eh9, ...
       integer_intervals_eh10, ...
       integer_intervals_eh11, ...
       ];


unique_intervals_eh = unique(sort(sum_intervals_eh));

Cov_eh11 = length(unique_intervals_eh) / 720 * 100;

%%
%% K_Satellite Setting(L0 / H12)

for n = 0:719
integer_intervals_kh1 = [];
integer_intervals_kh2 = [];
integer_intervals_kh3 = [];
integer_intervals_kh4 = [];
integer_intervals_kh5 = [];
integer_intervals_kh6 = [];
integer_intervals_kh7 = [];
integer_intervals_kh8 = [];
integer_intervals_kh9 = [];
integer_intervals_kh10 = [];
integer_intervals_kh11 = [];
integer_intervals_kh12 = [];

  
   for i = 1:length(ff_0)
   % 고궤도 satellite 1 위치 생성
       interval_kh1 = start_points1(i) :end_points1(i) ;
       interval_kh1(interval_kh1 > 719) = interval_kh1(interval_kh1 > 719) - 720; 
       integer_intervals_kh1 = [integer_intervals_kh1, interval_kh1];

   % 고궤도 satellite 2 위치 생성
       interval_kh2 = interval_kh1 + random_number_1; 
       interval_kh2(interval_kh2 > 719) = interval_kh2(interval_kh2 > 719) - 720;
       integer_intervals_kh2 = [integer_intervals_kh2, interval_kh2];

   % 고궤도 satellite 3 위치 생성
       interval_kh3 = interval_kh1 + random_number_2; 
       interval_kh3(interval_kh3 > 719) = interval_kh3(interval_kh3 > 719) - 720;
       integer_intervals_kh3 = [integer_intervals_kh3, interval_kh3];

   % 고궤도 satellite 4 위치 생성
       interval_kh4 = interval_kh1 + random_number_3; 
       interval_kh4(interval_kh4 > 719) = interval_kh4(interval_kh4 > 719) - 720;
       integer_intervals_kh4 = [integer_intervals_kh4, interval_kh4];

   % 고궤도 satellite 5 위치 생성
       interval_kh5 = interval_kh1 + random_number_4; 
       interval_kh5(interval_kh5 > 719) = interval_kh5(interval_kh5 > 719) - 720;
       integer_intervals_kh5 = [integer_intervals_kh5, interval_kh5];

   % 고궤도 satellite 6 위치 생성
       interval_kh6 = interval_kh1 + random_number_5; 
       interval_kh6(interval_kh6 > 719) = interval_kh6(interval_kh6 > 719) - 720;
       integer_intervals_kh6 = [integer_intervals_kh6, interval_kh6];

   % 고궤도 satellite 7 위치 생성
       interval_kh7 = interval_kh1 + random_number_6; 
       interval_kh7(interval_kh7 > 719) = interval_kh7(interval_kh7 > 719) - 720;
       integer_intervals_kh7 = [integer_intervals_kh7, interval_kh7];

   % 고궤도 satellite 8 위치 생성
       interval_kh8 = interval_kh1 + random_number_7; 
       interval_kh8(interval_kh8 > 719) = interval_kh8(interval_kh8 > 719) - 720;
       integer_intervals_kh8 = [integer_intervals_kh8, interval_kh8];

   % 고궤도 satellite 9 위치 생성
       interval_kh9 = interval_kh1 + random_number_8; 
       interval_kh9(interval_kh9 > 719) = interval_kh9(interval_kh9 > 719) - 720;
       integer_intervals_kh9 = [integer_intervals_kh9, interval_kh9];

   % 고궤도 satellite 10 위치 생성
       interval_kh10 = interval_kh1 + random_number_9; 
       interval_kh10(interval_kh10 > 719) = interval_kh10(interval_kh10 > 719) - 720;
       integer_intervals_kh10 = [integer_intervals_kh10, interval_kh10];

   % 고궤도 satellite 11 위치 생성
       interval_kh11 = interval_kh1 + random_number_10; 
       interval_kh11(interval_kh11 > 719) = interval_kh11(interval_kh11 > 719) - 720;
       integer_intervals_kh11 = [integer_intervals_kh11, interval_kh11];

   % 고궤도 satellite 12 랜덤 위치 생성
       interval_kh12 = interval_kh1 + n; 
       interval_kh12(interval_kh12 > 719) = interval_kh12(interval_kh12 > 719) - 720;
       integer_intervals_kh12 = [integer_intervals_kh12, interval_kh12];
   end


   sum_intervals_kh = ...
       [integer_intervals_kh1, ...
       integer_intervals_kh2, ...
       integer_intervals_kh3, ...
       integer_intervals_kh4, ...
       integer_intervals_kh5, ...
       integer_intervals_kh6, ...
       integer_intervals_kh7, ...
       integer_intervals_kh8, ...
       integer_intervals_kh9, ...
       integer_intervals_kh10, ...
       integer_intervals_kh11, ...
       integer_intervals_kh12, ...
       ];


unique_intervals_kh = unique(sort(sum_intervals_kh));

difference_kh = length(sum_intervals_kh) - length(unique_intervals_kh) ;

result_kh11 = [result_kh11; n, difference_kh];

end

%% E_Satellite Setting(L0 / H12)

for n = 0:719
integer_intervals_eh1 = [];
integer_intervals_eh2 = [];
integer_intervals_eh3 = [];
integer_intervals_eh4 = [];
integer_intervals_eh5 = [];
integer_intervals_eh6 = [];
integer_intervals_eh7 = [];
integer_intervals_eh8 = [];
integer_intervals_eh9 = [];
integer_intervals_eh10 = [];
integer_intervals_eh11 = [];
integer_intervals_eh12 = [];

  
   for i = 1:length(hh_0)
   % 고궤도 satellite 1 위치 생성
       interval_eh1 = start_points3(i) :end_points3(i) ;
       interval_eh1(interval_eh1 > 719) = interval_eh1(interval_eh1 > 719) - 720; 
       integer_intervals_eh1 = [integer_intervals_eh1, interval_eh1];

   % 고궤도 satellite 2 위치 생성
       interval_eh2 = interval_eh1 + random_number_1; 
       interval_eh2(interval_eh2 > 719) = interval_eh2(interval_eh2 > 719) - 720;
       integer_intervals_eh2 = [integer_intervals_eh2, interval_eh2];

   % 고궤도 satellite 3 위치 생성
       interval_eh3 = interval_eh1 + random_number_2; 
       interval_eh3(interval_eh3 > 719) = interval_eh3(interval_eh3 > 719) - 720;
       integer_intervals_eh3 = [integer_intervals_eh3, interval_eh3];

   % 고궤도 satellite 4 위치 생성
       interval_eh4 = interval_eh1 + random_number_3; 
       interval_eh4(interval_eh4 > 719) = interval_eh4(interval_eh4 > 719) - 720;
       integer_intervals_eh4 = [integer_intervals_eh4, interval_eh4];

   % 고궤도 satellite 5 위치 생성
       interval_eh5 = interval_eh1 + random_number_4; 
       interval_eh5(interval_eh5 > 719) = interval_eh5(interval_eh5 > 719) - 720;
       integer_intervals_eh5 = [integer_intervals_eh5, interval_eh5];

   % 고궤도 satellite 6 위치 생성
       interval_eh6 = interval_eh1 + random_number_5; 
       interval_eh6(interval_eh6 > 719) = interval_eh6(interval_eh6 > 719) - 720;
       integer_intervals_eh6 = [integer_intervals_eh6, interval_eh6];

   % 고궤도 satellite 7 위치 생성
       interval_eh7 = interval_eh1 + random_number_6; 
       interval_eh7(interval_eh7 > 719) = interval_eh7(interval_eh7 > 719) - 720;
       integer_intervals_eh7 = [integer_intervals_eh7, interval_eh7];

   % 고궤도 satellite 8 위치 생성
       interval_eh8 = interval_eh1 + random_number_7; 
       interval_eh8(interval_eh8 > 719) = interval_eh8(interval_eh8 > 719) - 720;
       integer_intervals_eh8 = [integer_intervals_eh8, interval_eh8];

   % 고궤도 satellite 9 위치 생성
       interval_eh9 = interval_eh1 + random_number_8; 
       interval_eh9(interval_eh9 > 719) = interval_eh9(interval_eh9 > 719) - 720;
       integer_intervals_eh9 = [integer_intervals_eh9, interval_eh9];

   % 고궤도 satellite 10 위치 생성
       interval_eh10 = interval_eh1 + random_number_9; 
       interval_eh10(interval_eh10 > 719) = interval_eh10(interval_eh10 > 719) - 720;
       integer_intervals_eh10 = [integer_intervals_eh10, interval_eh10];

   % 고궤도 satellite 11 위치 생성
       interval_eh11 = interval_eh1 + random_number_10; 
       interval_eh11(interval_eh11 > 719) = interval_eh11(interval_eh11 > 719) - 720;
       integer_intervals_eh11 = [integer_intervals_eh11, interval_eh11];

   % 고궤도 satellite 12 랜덤 위치 생성
       interval_eh12 = interval_eh1 + n; 
       interval_eh12(interval_eh12 > 719) = interval_eh12(interval_eh12 > 719) - 720;
       integer_intervals_eh12 = [integer_intervals_eh12, interval_eh12];
   end


   sum_intervals_eh = ...
       [integer_intervals_eh1, ...
       integer_intervals_eh2, ...
       integer_intervals_eh3, ...
       integer_intervals_eh4, ...
       integer_intervals_eh5, ...
       integer_intervals_eh6, ...
       integer_intervals_eh7, ...
       integer_intervals_eh8, ...
       integer_intervals_eh9, ...
       integer_intervals_eh10, ...
       integer_intervals_eh11, ...
       integer_intervals_eh12, ...
       ];


unique_intervals_eh = unique(sort(sum_intervals_eh));

difference_eh = length(sum_intervals_eh) - length(unique_intervals_eh) ;

result_eh11 = [result_eh11; n, difference_eh];

end

%% 공통

result_sum11 = [result_kh11(:,1), result_kh11(:,2) + result_eh11(:,2)];


rows_with_min_sum11 = result_sum11(result_sum11(:, 2) == min(result_sum11(:, 2)), :);

k_sum11 = rows_with_min_sum11(:,1);

random_index_11 = randi([1, length(k_sum11)]);

random_number_11 = k_sum11(random_index_11);

%% K_renewal 12
integer_intervals_kh1 = [];
integer_intervals_kh2 = [];
integer_intervals_kh3 = [];
integer_intervals_kh4 = [];
integer_intervals_kh5 = [];
integer_intervals_kh6 = [];
integer_intervals_kh7 = [];
integer_intervals_kh8 = [];
integer_intervals_kh9 = [];
integer_intervals_kh10 = [];
integer_intervals_kh11 = [];
integer_intervals_kh12 = [];

   for i = 1:length(ff_0)
   % 고궤도 satellite 1 위치 생성
       interval_kh1 = start_points1(i) :end_points1(i) ;
       interval_kh1(interval_kh1 > 719) = interval_kh1(interval_kh1 > 719) - 7320; 
       integer_intervals_kh1 = [integer_intervals_kh1, interval_kh1];

   % 고궤도 satellite 2 위치 재생성
       interval_kh2 = interval_kh1 + random_number_1; 
       interval_kh2(interval_kh2 > 719) = interval_kh2(interval_kh2 > 719) - 720;
       integer_intervals_kh2 = [integer_intervals_kh2, interval_kh2];

   % 고궤도 satellite 3 위치 재생성
       interval_kh3 = interval_kh1 + random_number_2; 
       interval_kh3(interval_kh3 > 719) = interval_kh3(interval_kh3 > 719) - 720;
       integer_intervals_kh3 = [integer_intervals_kh3, interval_kh3];

   % 고궤도 satellite 4 위치 재생성
       interval_kh4 = interval_kh1 + random_number_3; 
       interval_kh4(interval_kh4 > 719) = interval_kh4(interval_kh4 > 719) - 720;
       integer_intervals_kh4 = [integer_intervals_kh4, interval_kh4];

   % 고궤도 satellite 5 위치 재생성
       interval_kh5 = interval_kh1 + random_number_4; 
       interval_kh5(interval_kh5 > 719) = interval_kh5(interval_kh5 > 719) - 720;
       integer_intervals_kh5 = [integer_intervals_kh5, interval_kh5];

   % 고궤도 satellite 6 위치 재생성
       interval_kh6 = interval_kh1 + random_number_5; 
       interval_kh6(interval_kh6 > 719) = interval_kh6(interval_kh6 > 719) - 720;
       integer_intervals_kh6 = [integer_intervals_kh6, interval_kh6];

   % 고궤도 satellite 7 위치 재생성
       interval_kh7 = interval_kh1 + random_number_6; 
       interval_kh7(interval_kh7 > 719) = interval_kh7(interval_kh7 > 719) - 720;
       integer_intervals_kh7 = [integer_intervals_kh7, interval_kh7];

   % 고궤도 satellite 8 위치 재생성
       interval_kh8 = interval_kh1 + random_number_7; 
       interval_kh8(interval_kh8 > 719) = interval_kh8(interval_kh8 > 719) - 720;
       integer_intervals_kh8 = [integer_intervals_kh8, interval_kh8];

   % 고궤도 satellite 9 위치 재생성
       interval_kh9 = interval_kh1 + random_number_8; 
       interval_kh9(interval_kh9 > 719) = interval_kh9(interval_kh9 > 719) - 720;
       integer_intervals_kh9 = [integer_intervals_kh9, interval_kh9];

   % 고궤도 satellite 10 위치 재생성
       interval_kh10 = interval_kh1 + random_number_9; 
       interval_kh10(interval_kh10 > 719) = interval_kh10(interval_kh10 > 719) - 720;
       integer_intervals_kh10 = [integer_intervals_kh10, interval_kh10];

   % 고궤도 satellite 11 위치 재생성
       interval_kh11 = interval_kh1 + random_number_10; 
       interval_kh11(interval_kh11 > 719) = interval_kh11(interval_kh11 > 719) - 720;
       integer_intervals_kh11 = [integer_intervals_kh11, interval_kh11];

   % 고궤도 satellite 12 위치 재생성
       interval_kh12 = interval_kh1 + random_number_11; 
       interval_kh12(interval_kh12 > 719) = interval_kh12(interval_kh12 > 719) - 720;
       integer_intervals_kh12 = [integer_intervals_kh12, interval_kh12];
   end


   sum_intervals_kh = ...
       [integer_intervals_kh1, ...
       integer_intervals_kh2, ...
       integer_intervals_kh3, ...
       integer_intervals_kh4, ...
       integer_intervals_kh5, ...
       integer_intervals_kh6, ...
       integer_intervals_kh7, ...
       integer_intervals_kh8, ...
       integer_intervals_kh9, ...
       integer_intervals_kh10, ...
       integer_intervals_kh11, ...
       integer_intervals_kh12, ...
       ];


unique_intervals_kh = unique(sort(sum_intervals_kh));

Cov_kh12 = length(unique_intervals_kh) / 720 * 100;

%% E_renewal 11
integer_intervals_eh1 = [];
integer_intervals_eh2 = [];
integer_intervals_eh3 = [];
integer_intervals_eh4 = [];
integer_intervals_eh5 = [];
integer_intervals_eh6 = [];
integer_intervals_eh7 = [];
integer_intervals_eh8 = [];
integer_intervals_eh9 = [];
integer_intervals_eh10 = [];
integer_intervals_eh11 = [];
integer_intervals_eh12 = [];

   for i = 1:length(hh_0)
   % 고궤도 satellite 1 위치 생성
   interval_eh1 = start_points3(i) :end_points3(i) ;
       interval_eh1(interval_eh1 > 719) = interval_eh1(interval_eh1 > 719) - 720; 
       integer_intervals_eh1 = [integer_intervals_eh1, interval_eh1];

   % 고궤도 satellite 2 위치 재생성
       interval_eh2 = interval_eh1 + random_number_1; 
       interval_eh2(interval_eh2 > 719) = interval_eh2(interval_eh2 > 719) - 720;
       integer_intervals_eh2 = [integer_intervals_eh2, interval_eh2];

   % 고궤도 satellite 3 위치 재생성
       interval_eh3 = interval_eh1 + random_number_2; 
       interval_eh3(interval_eh3 > 719) = interval_eh3(interval_eh3 > 719) - 720;
       integer_intervals_eh3 = [integer_intervals_eh3, interval_eh3];

   % 고궤도 satellite 4 위치 재생성
       interval_eh4 = interval_eh1 + random_number_3; 
       interval_eh4(interval_eh4 > 719) = interval_eh4(interval_eh4 > 719) - 720;
       integer_intervals_eh4 = [integer_intervals_eh4, interval_eh4];

   % 고궤도 satellite 5 위치 재생성
       interval_eh5 = interval_eh1 + random_number_4; 
       interval_eh5(interval_eh5 > 719) = interval_eh5(interval_eh5 > 719) - 720;
       integer_intervals_eh5 = [integer_intervals_eh5, interval_eh5];

   % 고궤도 satellite 6 위치 재생성
       interval_eh6 = interval_eh1 + random_number_5; 
       interval_eh6(interval_eh6 > 719) = interval_eh6(interval_eh6 > 719) - 720;
       integer_intervals_eh6 = [integer_intervals_eh6, interval_eh6];

   % 고궤도 satellite 7 위치 재생성
       interval_eh7 = interval_eh1 + random_number_6; 
       interval_eh7(interval_eh7 > 719) = interval_eh7(interval_eh7 > 719) - 720;
       integer_intervals_eh7 = [integer_intervals_eh7, interval_eh7];

   % 고궤도 satellite 8 위치 재생성
       interval_eh8 = interval_eh1 + random_number_7; 
       interval_eh8(interval_eh8 > 719) = interval_eh8(interval_eh8 > 719) - 720;
       integer_intervals_eh8 = [integer_intervals_eh8, interval_eh8];

   % 고궤도 satellite 9 위치 재생성
       interval_eh9 = interval_eh1 + random_number_8; 
       interval_eh9(interval_eh9 > 719) = interval_eh9(interval_eh9 > 719) - 720;
       integer_intervals_eh9 = [integer_intervals_eh9, interval_eh9];

   % 고궤도 satellite 10 위치 재생성
       interval_eh10 = interval_eh1 + random_number_9; 
       interval_eh10(interval_eh10 > 719) = interval_eh10(interval_eh10 > 719) - 720;
       integer_intervals_eh10 = [integer_intervals_eh10, interval_eh10];

   % 고궤도 satellite 11 위치 재생성
       interval_eh11 = interval_eh1 + random_number_10; 
       interval_eh11(interval_eh11 > 719) = interval_eh11(interval_eh11 > 719) - 720;
       integer_intervals_eh11 = [integer_intervals_eh11, interval_eh11];

   % 고궤도 satellite 12 위치 재생성
       interval_eh12 = interval_eh1 + random_number_11; 
       interval_eh12(interval_eh12 > 719) = interval_eh12(interval_eh12 > 719) - 720;
       integer_intervals_eh12 = [integer_intervals_eh12, interval_eh12];
   end
   
   sum_intervals_eh = ...
       [integer_intervals_eh1, ...
       integer_intervals_eh2, ...
       integer_intervals_eh3, ...
       integer_intervals_eh4, ...
       integer_intervals_eh5, ...
       integer_intervals_eh6, ...
       integer_intervals_eh7, ...
       integer_intervals_eh8, ...
       integer_intervals_eh9, ...
       integer_intervals_eh10, ...
       integer_intervals_eh11, ...
       integer_intervals_eh12, ...
       ];


unique_intervals_eh = unique(sort(sum_intervals_eh));

Cov_eh12 = length(unique_intervals_eh) / 720 * 100;


%% 여기 전까지

random_N(1) = (0);
random_N(2) = (random_number_1);
random_N(3) = (random_number_2);
random_N(4) = (random_number_3);
random_N(5) = (random_number_4);
random_N(6) = (random_number_5);
random_N(7) = (random_number_6);
random_N(8) = (random_number_7);
random_N(9) = (random_number_8);
random_N(10) = (random_number_9);
random_N(11) = (random_number_10);
random_N(12) = (random_number_11)

percentage = 99.9;


Cov_kh12
Cov_eh12

  if Cov_kh12 >= percentage && Cov_eh12 >= percentage


        break;
    end
 trials = trials + 1; 
end

% 시뮬레이션 종료 시간 기록
elapsed_time = toc(start_time);

% 결과 출력
fprintf('시행 횟수: %d\n', trials);
fprintf('걸린 시간(초): %.2f\n', elapsed_time);



%%
%% 제궤도 그래프 K

n_1 = random_N +1;

%ff_0 = seedsatellite
ff_1 = ff_0 + n_1(1);
ff_2 = ff_0 + n_1(2);
ff_3 = ff_0 + n_1(3);
ff_4 = ff_0 + n_1(4);
ff_5 = ff_0 + n_1(5);
ff_6 = ff_0 + n_1(6);
ff_7 = ff_0 + n_1(7);
ff_8 = ff_0 + n_1(8);
ff_9 = ff_0 + n_1(9);
ff_10 = ff_0 + n_1(10);
ff_11 = ff_0 + n_1(11);
ff_12 = ff_0 + n_1(12);


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




%% 고궤도 그래프 K

n_2 = [0] +1;

%gg_0 = seedsatellite
gg_1 = gg_0 + n_2(1);
% gg_2 = gg_0 + n_2(2);
% gg_3 = gg_0 + n_2(3);
% gg_4 = gg_0 + n_2(4);
% gg_5 = gg_0 + n_2(5);
% gg_6 = gg_0 + n_2(6);
% gg_7 = gg_0 + n_2(7);


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
RV2 = [];
current_value = all_x_coords2(1);
start_index = 1;

% 각 원소에 대해 반복하면서 연속된 값의 구간을 찾음
for i = 2:length(all_x_coords2)
    if all_x_coords2(i) ~= current_value
        % 현재 값과 다른 값이 나타난 경우, 해당 구간의 시작과 끝을 저장하고 현재 값 갱신
        end_index = i - 1;
        RV2 = [RV2; current_value, end_index - start_index + 1];
        current_value = all_x_coords2(i);
        start_index = i;
    end
end

% 마지막 구간에 대한 처리
end_index = length(all_x_coords2);
RV2 = [RV2; current_value, end_index - start_index+1];


%% 빠진 숫자 찾기

% 주어진 행렬
matrix2 = RV2 ;

% 현재 1열에서 빠진 숫자를 찾기
missing_numbers = setdiff(0:L-1, matrix2(:,1));

% 새롭게 추가된 행에 대해 0으로 입력
new_rows = zeros(length(missing_numbers), 2);
new_rows(:, 1) = missing_numbers;

% 새로운 행렬 생성
new_matrix2 = [matrix2; new_rows];


% 1열을 기준으로 행렬을 정렬
sorted_matrix2 = sortrows(new_matrix2, 1);

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


%% K

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
sorted_matrix3(:, 1) = sorted_matrix3(:, 1);  % 1열은 첫 번째 행렬에서 복사

% 2열의 OR 연산을 사용하여 값 결정
sorted_matrix3(:, 2) = sorted_matrix3(:, 2) | sorted_matrix3(:, 2);

sorted_matrix3(:,2) = sorted_matrix3(:,2) + sorted_matrix3(:,2);

count1 = sum(sorted_matrix1(:, 2) >= 1); 
count2 = sum(sorted_matrix2(:, 2) >= 1); 
count = sum(sorted_matrix3(:, 2) >= 1); 

coverage_total = (count) / L * 100;

all_x_coords3 = unique(sort(all_x_coords1(:,2) + all_x_coords2(:,2)));

%%
%% 제궤도 그래프 E

n_3 = random_N +1;

%hh_0 = seedsatellite
hh_1 = hh_0 + n_3(1);
hh_2 = hh_0 + n_3(2);
hh_3 = hh_0 + n_3(3);
hh_4 = hh_0 + n_3(4);
hh_5 = hh_0 + n_3(5);
hh_6 = hh_0 + n_3(6);
hh_7 = hh_0 + n_3(7);
hh_8 = hh_0 + n_3(8);
hh_9 = hh_0 + n_3(9);
hh_10 = hh_0 + n_3(10);
hh_11 = hh_0 + n_3(11);
hh_12 = hh_0 + n_3(12);

% x축 범위 설정
x_range = 0:L-1;

for idx = 1:length(n_3)
    hh = eval(['hh_', num2str(idx)]);
    x_coords3{idx} = []; % 해당 작업공간의 x좌표를 저장할 배열 초기화

    for i = 1:size(hh,1)
        start_time = hh(i, 1);
        end_time = hh(i, 2);

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
          x_coords3{idx} = [x_coords3{idx}, x_range(start_index:L),(0:end_index)];

        else
          x_coords3{idx} = [x_coords3{idx}, x_range(start_index:end_index)];
        end



    end               
end

% 모든 작업공간의 x 좌표를 모은 배열 초기화
all_x_coords3 = [];

%% 각 작업공간의 x 좌표를 all_x_coords3에 추가하고 작은 수부터 오름차순으로 정렬
for idx = 1:length(n_3)
    all_x_coords3 = [all_x_coords3, x_coords3{idx}];
end

all_x_coords3 = sort(all_x_coords3);

% 중복된 값을 제거하고 유일한 값들만 추출
unique_x_coords3 = unique(all_x_coords3);

coverage4 = length(unique_x_coords3) / L * 100;

%%

% 초기화
RV3 = [];
current_value = all_x_coords3(1);
start_index = 1;

% 각 원소에 대해 반복하면서 연속된 값의 구간을 찾음
for i = 2:length(all_x_coords3)
    if all_x_coords3(i) ~= current_value
        % 현재 값과 다른 값이 나타난 경우, 해당 구간의 시작과 끝을 저장하고 현재 값 갱신
        end_index = i - 1;
        RV3 = [RV3; current_value, end_index - start_index + 1];
        current_value = all_x_coords3(i);
        start_index = i;
    end
end

% 마지막 구간에 대한 처리
end_index = length(all_x_coords3);
RV3 = [RV3; current_value, end_index - start_index+1];


%% 빠진 숫자 찾기

% 주어진 행렬
matrix3 = RV3 ;

% 현재 1열에서 빠진 숫자를 찾기
missing_numbers = setdiff(0:L-1, matrix3(:,1));

% 새롭게 추가된 행에 대해 0으로 입력
new_rows = zeros(length(missing_numbers), 2);
new_rows(:, 1) = missing_numbers;

% 새로운 행렬 생성
new_matrix3 = [matrix3; new_rows];


% 1열을 기준으로 행렬을 정렬
sorted_matrix3 = sortrows(new_matrix3, 1);

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
for i = 1:size(ii_0, 1)
    if neg3(i) < 0
        start_time_new = L + neg3(i,1);
        end_time_new = L + neg3(i,2);
    else 
        start_time_new = hh_0(i, 1);
        end_time_new = hh_0(i, 2);
    end
    start_time = hh_0(i, 1);
    end_time = hh_0(i, 2);

    if end_time_new > L-1
        end_time_new = L-1;
    end

    % 시작 시간과 종료 시간에 해당하는 x 인덱스 찾기
    start_index = find(x_range >= start_time, 1);
    end_index = find(x_range <= end_time, 1, 'last');
    stn = start_time_new;
    etn = end_time_new;

    % y 값 설정
    y_values3a = zeros(size(x_range));
    y_values3a(start_index:end_index) = 1;
    if neg3(i) < 0
        y_values3a(stn:etn) = 1;
    end

    % 그래프 그리기
    plot(x_range, y_values3a, '-o r', 'LineWidth', 1, 'MarkerSize', 1, 'Color', [red1, green1, blue1]);
    hold on;
    area(x_range, y_values3a, 'FaceColor', [red11, green11, blue11]);
    plot(x_range, y_values3a, 'Color', [red1, green1, blue1], 'LineWidth', 1.5); 
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
    y_values3p = zeros(size(x_range));
    y_values3p(n_3(1,:)) = 1;


plot(x_range, y_values3p, 'Color', [red1, green1, blue1], 'LineWidth', 1, 'MarkerSize', 1);
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
plot(sorted_matrix3(:,1), sorted_matrix3(:,2), 'Color', [red1, green1, blue1], 'LineWidth', 1, 'MarkerSize', 1);
title('Coverage Timeline');
hold on;
area(sorted_matrix3(:,1), sorted_matrix3(:,2), 'FaceColor', [red11, green11, blue11]);
plot(sorted_matrix3(:,1), sorted_matrix3(:,2), 'Color', [red1, green1, blue1], 'LineWidth', 1);

xlim([0, L-1]);
ylim([0, 3]);

% 그래프 세부 설정
xlabel('n');
ylabel('bj[n]');
grid on;
xlim([0, L-1]);




%% 고궤도 그래프 E

n_4 = [0] +1;

%hh_0 = seedsatellite
ii_1 = ii_0 + n_4(1);
% ii_2 = ii_0 + n_4(2);
% ii_3 = ii_0 + n_4(3);
% ii_4 = ii_0 + n_4(4);
% ii_5 = ii_0 + n_4(5);
% ii_6 = ii_0 + n_4(6);
% ii_7 = ii_0 + n_4(7);


% x축 범위 설정
x_range = 0:L-1;

for idx = 1:length(n_4)
    ii = eval(['ii_', num2str(idx)]);
    x_coords4{idx} = []; % 해당 작업공간의 x좌표를 저장할 배열 초기화

    for i = 1:size(ii,1)
        start_time = ii(i, 1);
        end_time = ii(i, 2);

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
          x_coords4{idx} = [x_coords4{idx}, x_range(start_index:L),(0:end_index)];

        else
          x_coords4{idx} = [x_coords4{idx}, x_range(start_index:end_index)];
        end



    end               
end

% 모든 작업공간의 x 좌표를 모은 배열 초기화
all_x_coords4 = [];

%% 각 작업공간의 x 좌표를 all_x_coords3에 추가하고 작은 수부터 오름차순으로 정렬
for idx = 1:length(n_4)
    all_x_coords4 = [all_x_coords4, x_coords4{idx}];
end

all_x_coords4 = sort(all_x_coords4);

% 중복된 값을 제거하고 유일한 값들만 추출
unique_x_coords4 = unique(all_x_coords4);

coverage4 = length(unique_x_coords4) / L * 100;

%%

% 초기화
RV4 = [];
current_value = all_x_coords4(1);
start_index = 1;

% 각 원소에 대해 반복하면서 연속된 값의 구간을 찾음
for i = 2:length(all_x_coords4)
    if all_x_coords4(i) ~= current_value
        % 현재 값과 다른 값이 나타난 경우, 해당 구간의 시작과 끝을 저장하고 현재 값 갱신
        end_index = i - 1;
        RV4 = [RV4; current_value, end_index - start_index + 1];
        current_value = all_x_coords4(i);
        start_index = i;
    end
end

% 마지막 구간에 대한 처리
end_index = length(all_x_coords4);
RV4 = [RV4; current_value, end_index - start_index+1];


%% 빠진 숫자 찾기

% 주어진 행렬
matrix4 = RV4 ;

% 현재 1열에서 빠진 숫자를 찾기
missing_numbers = setdiff(0:L-1, matrix4(:,1));

% 새롭게 추가된 행에 대해 0으로 입력
new_rows = zeros(length(missing_numbers), 2);
new_rows(:, 1) = missing_numbers;

% 새로운 행렬 생성
new_matrix4 = [matrix4; new_rows];


% 1열을 기준으로 행렬을 정렬
sorted_matrix4 = sortrows(new_matrix4, 1);

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
for i = 1:size(hh_0, 1)
    if neg4(i) < 0
        start_time_new = L + neg4(i,1);
        end_time_new = L + neg4(i,2);
    else 
        start_time_new = ii_0(i, 1);
        end_time_new = ii_0(i, 2);
    end
    start_time = ii_0(i, 1);
    end_time = ii_0(i, 2);

    if end_time_new > L-1
        end_time_new = L-1;
    end

    % 시작 시간과 종료 시간에 해당하는 x 인덱스 찾기
    start_index = find(x_range >= start_time, 1);
    end_index = find(x_range <= end_time, 1, 'last');
    stn = start_time_new;
    etn = end_time_new;

    % y 값 설정
    y_values4a = zeros(size(x_range));
    y_values4a(start_index:end_index) = 1;
    if neg4(i) < 0
        y_values4a(stn:etn) = 1;
    end

    % 그래프 그리기
    plot(x_range, y_values4a, '-o r', 'LineWidth', 1, 'MarkerSize', 1, 'Color', [red2, green2, blue2]);
    hold on;
    area(x_range, y_values4a, 'FaceColor', [red22, green22, blue22]);
    plot(x_range, y_values4a, 'Color', [red2, green2, blue2], 'LineWidth', 1.5); 
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
    y_values4p = zeros(size(x_range));
    y_values4p(n_4(1,:)) = 1;


plot(x_range, y_values4p, 'Color', [red2, green2, blue2], 'LineWidth', 1, 'MarkerSize', 1);
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
plot(sorted_matrix4(:,1), sorted_matrix4(:,2), 'Color', [red2, green2, blue2], 'LineWidth', 1, 'MarkerSize', 1);
title('Coverage Timeline');
hold on;
area(sorted_matrix4(:,1), sorted_matrix4(:,2), 'FaceColor', [red22, green22, blue22]);
plot(sorted_matrix4(:,1), sorted_matrix4(:,2), 'Color', [red2, green2, blue2], 'LineWidth', 1);

xlim([0, L-1]);
ylim([0, 3]);

% 그래프 세부 설정
xlabel('n');
ylabel('bj[n]');
grid on;
xlim([0, L-1]);


%% E

% 새로운 그래프에 두 개의 플롯을 함께 보기
figure3 = figure;
subplot(3, 1, 1);

hold on;


% hh
for i = 1:size(hh_0, 1)
    if neg3(i) < 0
        start_time_new = L + neg3(i,1);
        end_time_new = L + neg3(i,2);
    else 
        start_time_new = hh_0(i, 1);
        end_time_new = hh_0(i, 2);
    end
    start_time = hh_0(i, 1);
    end_time = hh_0(i, 2);

    if end_time_new > L-1
        end_time_new = L-1;
    end

    % 시작 시간과 종료 시간에 해당하는 x 인덱스 찾기
    start_index = find(x_range >= start_time, 1);
    end_index = find(x_range <= end_time, 1, 'last');
    stn = start_time_new;
    etn = end_time_new;

    % y 값 설정
    y_values3a = zeros(size(x_range));
    y_values3a(start_index:end_index) = 1;
    if neg3(i) < 0
        y_values3a(stn:etn) = 1;
    end

    % 그래프 그리기
    hold on;
    plot(x_range, y_values3a, '-o r', 'LineWidth', 1, 'MarkerSize', 1, 'Color', [red1, green1, blue1]);
    area(x_range, y_values3a, 'FaceColor', [red11, green11, blue11], 'FaceAlpha', 0.7);
    plot(x_range, y_values3a, 'Color', [red1, green1, blue1], 'LineWidth', 1.5); 
end

% ii
for i = 1:size(ii_0, 1)
    if neg4(i) < 0
        start_time_new = L + neg4(i,1);
        end_time_new = L + neg4(i,2);
    else 
        start_time_new = ii_0(i, 1);
        end_time_new = ii_0(i, 1);
    end
    start_time = ii_0(i, 1);
    end_time = ii_0(i, 2);

    if end_time_new > L-1
        end_time_new = L-1;
    end

    % 시작 시간과 종료 시간에 해당하는 x 인덱스 찾기
    start_index = find(x_range >= start_time, 1);
    end_index = find(x_range <= end_time, 1, 'last');
    stn = start_time_new;
    etn = end_time_new;

    % y 값 설정
    y_values4a = zeros(size(x_range));
    y_values4a(start_index:end_index) = 1;
    if neg4(i) < 0
        y_values4a(stn:etn) = 1;
    end

    % 그래프 그리기
    plot(x_range, y_values4a, '-o r', 'LineWidth', 1, 'MarkerSize', 1, 'Color', [red2, green2, blue2]);
    hold on;
    area(x_range, y_values4a, 'FaceColor', [red22, green22, blue22], 'FaceAlpha', 0.7);
    plot(x_range, y_values4a, 'Color', [red2, green2, blue2], 'LineWidth', 1.5);
end


hold on; % 여러 플롯을 함께 그리기 위해 hold on 사용
plot(x_range, y_values3a, 'Color', [red1, green1, blue1]); % 첫 번째 플롯 복사
plot(x_range, y_values4a, 'Color', [red2, green2, blue2]); % 두 번째 플롯 복사
xlim([0, L-1]);
ylim([0, 1]);
yticks([0, 1]);

% 그래프 세부 설정
ylabel('vo,j[n]');
title('Seed Satellite Access Profile');
grid on;



subplot(3, 1, 2);
hold on; % 여러 플롯을 함께 그리기 위해 hold on 사용
plot(x_range, y_values3p, 'Color', [red1, green1, blue1]); % 첫 번째 플롯 복사
plot(x_range, y_values4p, 'Color', [red2, green2, blue2]); % 두 번째 플롯 복사
xlim([0, L-1]);
ylim([0, 1]);
yticks([0, 1]);

% 그래프 세부 설정
ylabel('x[n]');
title('Constellation Pattern Vector');
grid on;



subplot(3, 1, 3);
hold on; % 여러 플롯을 함께 그리기 위해 hold on 사용
plot(sorted_matrix3(:,1), sorted_matrix3(:,2), 'Color', [red1, green1, blue1]); % 첫 번째 플롯 복사
hold on;
area(sorted_matrix3(:,1), sorted_matrix3(:,2), 'FaceColor', [red11, green11, blue11], 'FaceAlpha', 0.7);
plot(sorted_matrix3(:,1), sorted_matrix3(:,2), 'Color', [red1, green1, blue1], 'LineWidth', 1.5);

plot(sorted_matrix4(:,1), sorted_matrix4(:,2), 'Color', [red2, green2, blue2]); % 두 번째 플롯 복사
hold on;
area(sorted_matrix4(:,1), sorted_matrix4(:,2), 'FaceColor', [red22, green22, blue22], 'FaceAlpha', 0.7);
plot(sorted_matrix4(:,1), sorted_matrix4(:,2), 'Color', [red2, green2, blue2], 'LineWidth', 1.5);

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
sorted_matrix3(:, 1) = sorted_matrix3(:, 1);  % 1열은 첫 번째 행렬에서 복사

% 2열의 OR 연산을 사용하여 값 결정
sorted_matrix3(:, 2) = sorted_matrix3(:, 2) | sorted_matrix4(:, 2);

sorted_matrix3(:,2) = sorted_matrix3(:,2) + sorted_matrix4(:,2);

count1 = sum(sorted_matrix3(:, 2) >= 1); 
count2 = sum(sorted_matrix4(:, 2) >= 1); 
count = sum(sorted_matrix3(:, 2) >= 1); 

coverage_total = (count) / L * 100;

all_x_coords3 = unique(sort(all_x_coords3(:,2) + all_x_coords4(:,2)));
