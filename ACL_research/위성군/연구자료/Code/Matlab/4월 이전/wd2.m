function walkerDeltaScenario()

    % Walker Delta 배열 파라미터 설정
    N = 3; % 각 궤도 평면에 있는 위성 수
    M = 3; % 총 궤도 평면 수
    semiMajorAxis = 7000e3; % 궤도의 반장축 (예시 값)
    inclination = 45; % 궤도의 경사각 (예시 값)

    % 시나리오 생성
    scenario = aero.Scenario('Name', 'Walker Delta Scenario');

    % Walker Delta 배열에 따라 위성 생성 및 시나리오에 추가
    for m = 1:M
        for n = 1:N
            % 각 위성의 초기 위상을 설정 (등간격으로 설정)
            trueAnomaly = 360 / N * (n - 1);

            % Walker Delta 배열에 따른 위성 생성 및 추가
            sat = satellite(scenario, semiMajorAxis, 0, inclination, 0, 0, trueAnomaly);
        end
    end

    % 시나리오 뷰어 열기
    aero.view(scenario);

end
