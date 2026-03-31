clc
clear all


% 주어진 행렬 ff_0
ff_0 = [292	294
340	345
391	395
442	445
493	496
542	547
594	595];



% 구간의 시작점과 끝점을 변수에 저장
start_points = ff_0(:, 1);
end_points = ff_0(:, 2);

% 모든 구간을 정수로 표시해서 작업 공간에 저장
result = [];

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
integer_intervals_20 = [];
integer_intervals_21 = [];
integer_intervals_22 = [];
integer_intervals_23 = [];
integer_intervals_24 = [];
integer_intervals_25 = [];
integer_intervals_26 = [];
integer_intervals_27 = [];
integer_intervals_28 = [];
integer_intervals_29 = [];
integer_intervals_30 = [];
integer_intervals_31 = [];
integer_intervals_32 = [];
integer_intervals_33 = [];
integer_intervals_34 = [];

%%

for i = 1:length(ff_0)
    interval_1 = start_points(i):end_points(i) ;
    interval_1(interval_1 > 719) = interval_1(interval_1 > 719) - 720; 
    integer_intervals_1 = [integer_intervals_1, interval_1];

    interval_2 = interval_1 + 6; 
    interval_2(interval_2 > 719) = interval_2(interval_2 > 719) - 720;
    integer_intervals_2 = [integer_intervals_2, interval_2];

    interval_3 = interval_1 + 12; 
    interval_3(interval_3 > 719) = interval_3(interval_3 > 719) - 720;
    integer_intervals_3 = [integer_intervals_3, interval_3];

    interval_4 = interval_1 + 18; 
    interval_4(interval_4 > 719) = interval_4(interval_4 > 719) - 720;
    integer_intervals_4 = [integer_intervals_4, interval_4];

    interval_5 = interval_1 + 24; 
    interval_5(interval_5 > 719) = interval_5(interval_5 > 719) - 720;
    integer_intervals_5 = [integer_intervals_5, interval_5];

    interval_6 = interval_1 + 30; 
    interval_6(interval_6 > 719) = interval_6(interval_6 > 719) - 720;
    integer_intervals_6 = [integer_intervals_6, interval_6];

    interval_7 = interval_1 + 36; 
    interval_7(interval_7 > 719) = interval_7(interval_7 > 719) - 720;
    integer_intervals_7 = [integer_intervals_7, interval_7];

    interval_8 = interval_1 + 42; 
    interval_8(interval_8 > 719) = interval_8(interval_8 > 719) - 720;
    integer_intervals_8 = [integer_intervals_8, interval_8];

    interval_9 = interval_1 + 298; 
    interval_9(interval_9 > 719) = interval_9(interval_9 > 719) - 720;
    integer_intervals_9 = [integer_intervals_9, interval_9];

    interval_10 = interval_1 + 304; 
    interval_10(interval_10 > 719) = interval_10(interval_10 > 719) - 720;
    integer_intervals_10 = [integer_intervals_10, interval_10];

    interval_11 = interval_1 + 310; 
    interval_11(interval_11 > 719) = interval_11(interval_11 > 719) - 720;
    integer_intervals_11 = [integer_intervals_11, interval_11];

    interval_12 = interval_1 + 316; 
    interval_12(interval_12 > 719) = interval_12(interval_12 > 719) - 720;
    integer_intervals_12 = [integer_intervals_12, interval_12];

    interval_13 = interval_1 + 322; 
    interval_13(interval_13 > 719) = interval_13(interval_13 > 719) - 720;
    integer_intervals_13 = [integer_intervals_13, interval_13];

    interval_14 = interval_1 + 328; 
    interval_14(interval_14 > 719) = interval_14(interval_14 > 719) - 720;
    integer_intervals_14 = [integer_intervals_14, interval_14];

    interval_15 = interval_1 + 334; 
    interval_15(interval_15 > 719) = interval_15(interval_15 > 719) - 720;
    integer_intervals_15 = [integer_intervals_15, interval_15];

    interval_16 = interval_1 + 340; 
    interval_16(interval_16 > 719) = interval_16(interval_16 > 719) - 720;
    integer_intervals_16 = [integer_intervals_16, interval_16];

    interval_17 = interval_1 + 395; 
    interval_17(interval_17 > 719) = interval_17(interval_17 > 719) - 720;
    integer_intervals_17 = [integer_intervals_17, interval_17];

    interval_18 = interval_1 + 614; 
    interval_18(interval_18 > 719) = interval_18(interval_18 > 719) - 720;
    integer_intervals_18 = [integer_intervals_18, interval_18];

    interval_19 = interval_1 + 452; 
    interval_19(interval_19 > 719) = interval_19(interval_19 > 719) - 720;
    integer_intervals_19 = [integer_intervals_19, interval_19];

    interval_20 = interval_1 + 457; 
    interval_20(interval_20 > 719) = interval_20(interval_20 > 719) - 720;
    integer_intervals_20 = [integer_intervals_20, interval_20];

    interval_21 = interval_1 + 419; 
    interval_21(interval_21 > 719) = interval_21(interval_21 > 719) - 720;
    integer_intervals_21 = [integer_intervals_21, interval_21];

    interval_22 = interval_1 + 425; 
    interval_22(interval_22 > 719) = interval_22(interval_22 > 719) - 720;
    integer_intervals_22 = [integer_intervals_22, interval_22];

    interval_23 = interval_1 + 431; 
    interval_23(interval_23 > 719) = interval_23(interval_23 > 719) - 720;
    integer_intervals_23 = [integer_intervals_23, interval_23];

    interval_24 = interval_1 + 437; 
    interval_24(interval_24 > 719) = interval_24(interval_24 > 719) - 720;
    integer_intervals_24 = [integer_intervals_24, interval_24];

    interval_25 = interval_1 + 647; 
    interval_25(interval_25 > 719) = interval_25(interval_25 > 719) - 720;
    integer_intervals_25 = [integer_intervals_25, interval_25];

    interval_26 = interval_1 + 641; 
    interval_26(interval_26 > 719) = interval_26(interval_26 > 719) - 720;
    integer_intervals_26 = [integer_intervals_26, interval_26];

    interval_27 = interval_1 + 3; 
    interval_27(interval_27 > 719) = interval_27(interval_27 > 719) - 720;
    integer_intervals_27 = [integer_intervals_27, interval_27];

    interval_28 = interval_1 + 172; 
    interval_28(interval_28 > 719) = interval_28(interval_28 > 719) - 720;
    integer_intervals_28 = [integer_intervals_28, interval_28];

    interval_29 = interval_1 + 9; 
    interval_29(interval_29 > 719) = interval_29(interval_29 > 719) - 720;
    integer_intervals_29 = [integer_intervals_29, interval_29];

    interval_30 = interval_1 + 15; 
    interval_30(interval_30 > 719) = interval_30(interval_30 > 719) - 720;
    integer_intervals_30 = [integer_intervals_30, interval_30];

    interval_31 = interval_1 + 148; 
    interval_31(interval_31 > 719) = interval_31(interval_31 > 719) - 720;
    integer_intervals_31 = [integer_intervals_31, interval_31];

    interval_32 = interval_1 + 416; 
    interval_32(interval_32 > 719) = interval_32(interval_32 > 719) - 720;
    integer_intervals_32 = [integer_intervals_32, interval_32];

    interval_33 = interval_1 + 40; 
    interval_33(interval_33 > 719) = interval_33(interval_33 > 719) - 720;
    integer_intervals_33 = [integer_intervals_33, interval_33];

    interval_34 = interval_1 + n; 
    interval_34(interval_34 > 719) = interval_34(interval_34 > 719) - 720;
    integer_intervals_34 = [integer_intervals_34, interval_34];
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
    integer_intervals_23, ...
    integer_intervals_24, ...
    integer_intervals_25, ...
    integer_intervals_26, ...
    integer_intervals_27, ...
    integer_intervals_28, ...
    integer_intervals_29, ...
    integer_intervals_30, ...
    integer_intervals_31, ...
    integer_intervals_32, ...
    integer_intervals_33, ...
    integer_intervals_34, ...
    ];

% 행렬을 오름차순으로 정렬하고 중복된 요소를 제거합니다.
unique_intervals = unique(sort(sum_intervals));


% 두 행렬의 길이 차이 계산
difference = length(sum_intervals) - length(unique_intervals) ;

 % 결과를 작업 공간에 저장
    result = [result; n, difference];
end


A = length(unique_intervals) / 720 * 100

% result 변수에서 2열이 최소인 행을 찾습니다.
rows_with_min = result(result(:, 2) == min(result(:, 2)), :);

k_1 = rows_with_min(:,1);
