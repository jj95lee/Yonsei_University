% 누적된 행렬을 반환하는 함수
function accumulator = accumulate_results()
    % 1x10 랜덤 행렬 생성
    new_row = randi(100, 1, 10); % 1부터 100 사이의 랜덤 정수 10개
    
    % 지속적인 변수를 선언
    persistent accumulated_data;
    
    % 만약 첫 번째 실행이라면, 새로운 행렬 생성
    if isempty(accumulated_data)
        accumulated_data = new_row; % 초기화
    else
        % 기존 데이터에 새 행 추가
        accumulated_data = [accumulated_data; new_row];
    end
    
    % 결과 반환
    accumulator = accumulated_data;
end
