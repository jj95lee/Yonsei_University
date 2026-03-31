% 가능한 모든 조합 생성
combinations = dec2bin(0:2^9-1) - '0';
selected_combinations = combinations(sum(combinations, 2) >= 1, :);

% 각 조합의 합 계산
sums = sum(selected_combinations, 2);

% 모든 성분이 1 이상인 조합 선택
valid_combinations = selected_combinations(sums >= 1, :);

disp('모든 성분이 1 이상인 최소 조합:');
disp(valid_combinations);
disp(['조합의 수: ', num2str(size(valid_combinations, 1))]);
