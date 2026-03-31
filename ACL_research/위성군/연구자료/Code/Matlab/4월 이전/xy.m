% 예시 데이터 생성 (실제 데이터로 대체해야 함)
satellite_data = ab;

% 중복 제거
unique_satellite_data = unique(satellite_data, 'rows');

% 그래프 초기화
figure;

% 그래프에 연결된 시간 표시
for i = 1:size(unique_satellite_data, 1)
    start_time = unique_satellite_data(i, 1);
    end_time = unique_satellite_data(i, 2);
    
    % 연결된 구간을 1로 설정
    plot([start_time, end_time], [1, 1], 'LineWidth', 2, 'Color', 'b');
    hold on;

    
end

% 그래프 스타일 및 레이블 설정
xlabel('Time (seconds)');
ylabel('Satellite Connection');
title('Connected Time Intervals');

% x축 범위 설정
xlim([0, 86400]);