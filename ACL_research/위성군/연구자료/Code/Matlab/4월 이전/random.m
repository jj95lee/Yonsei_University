% 행렬 크기 정의
row_count = 72;
col_count = 1;

% 목표 1의 비율
target_ratio = 0.78;

% 초기 행렬 생성
Y = randi([0, 1], row_count, col_count);

% 초기 1의 비율 계산
current_ratio = sum(Y) / row_count;

% 1의 비율이 목표에 미치지 못할 때까지 반복
iteration_count = 0;
while current_ratio < target_ratio
    % 새로운 행렬 생성
    Y = randi([0, 1], row_count, col_count);
    
    % 1의 비율 업데이트
    current_ratio = sum(Y) / row_count;
    
    % 시행 횟수 증가
    iteration_count = iteration_count + 1;
end

% 결과 출력
disp(['1의 비율이 목표인 ', num2str(target_ratio * 100), '% 이상이 되기까지 ', num2str(iteration_count), '번 시도했습니다.']);
