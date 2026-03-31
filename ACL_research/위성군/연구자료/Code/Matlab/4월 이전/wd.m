clear all


startTime = datetime('2024-01-01 00:00:00');
stopTime = datetime('2024-01-02 00:00:00');
sampleTime = 120; % seconds
sc = satelliteScenario(startTime,stopTime,sampleTime)

sat(1) = walkerDelta(sc,7270000, 45, 6, 2, 1, "ArgumentOfLatitude", 0, "OrbitPropagator","two-body-keplerian")

names = sat.Name + " Camera";
cam = conicalSensor(sat,"Name",names,"MaxViewAngle",179.9999)

name = "Ground Station";
latitude = 43.1; % 위도
longitude = 131.9; % 경도
altitude = 0; % 고도
minElevationAngle = 0; % degrees
geoSite = groundStation(sc, ...
    "Name",name,"Latitude", latitude, "Longitude", longitude, "Altitude", altitude, ...
     "MinElevationAngle",minElevationAngle);



ac = access(sat,geoSite);

% Properties of access analysis objects
ac(1);

v = satelliteScenarioViewer(sc,"ShowDetails",true);
sat.ShowLabel = true;
geoSite.ShowLabel = true;
show(sat);

%fov = fieldOfView(cam)
% fov = fieldOfView(cam([cam.Name] == "Satellite 2 Camera"))
% fov = fieldOfView(cam([cam.Name] == "Satellite 3 Camera"))
% fov = fieldOfView(cam([cam.Name] == "Satellite 4 Camera"))
% fov = fieldOfView(cam([cam.Name] == "Satellite 5 Camera"))
% fov = fieldOfView(cam([cam.Name] == "Satellite 6 Camera"))
% fov = fieldOfView(cam([cam.Name] == "Satellite 7 Camera"))
% fov = fieldOfView(cam([cam.Name] == "Satellite 8 Camera"))
% fov = fieldOfView(cam([cam.Name] == "Satellite 9 Camera"))
% fov = fieldOfView(cam([cam.Name] == "Satellite 10 Camera"))

accessIntervals(ac);

play(sc);

for idx = 1:numel(ac)
    [s,time] = accessStatus(ac(idx));
    
    if idx == 1
        % Initialize system-wide access status vector in the first iteration
        systemWideAccessStatus = s;
    else
        % Update system-wide access status vector by performing a logical OR
        % with access status for the current camera-site access
        % analysis
        systemWideAccessStatus = or(systemWideAccessStatus,s);
    end
end

% 총 시간 (초)과 총 x 개수 설정
totalSeconds = seconds(sc.StopTime - sc.StartTime);
x_count = 720;

% 0초부터 120초 간격으로 720개의 x값 생성
x_seconds = 0:sampleTime:(x_count - 1) * sampleTime;

% systemWideAccessStatus 벡터를 720개로 줄이기
systemWideAccessStatus = systemWideAccessStatus(1:720);

% 그래프 다시 그리기
plot(x_seconds/120, systemWideAccessStatus, "LineWidth", 2);
xlim([0 720]);
grid on;
xlabel("n");  % 초 단위로 x축 표시
ylabel("v0j[n]");
title('Seed satellite access profile')

% 아래 그래프를 색칠할 때 사용할 벡터 생성
x_seconds_shaded = x_seconds/120; % x축 벡터

% System-Wide Access Status 벡터를 0 또는 1 값으로 변환
shaded_values = systemWideAccessStatus; % 여기에 색칠할 영역의 값이 있어야 함

% 1인 값(접근 가능한 상태)을 0.5로 변환하여 색칠
shaded_values(shaded_values == 1) = 0.5;

% 그래프 그리기
plot(x_seconds_shaded, systemWideAccessStatus, "LineWidth", 2);
xlim([0 x_count]);
grid on;
xlabel("n");  % 초 단위로 x축 표시
ylabel("v0j[n]");
title('Seed satellite access profile')

% 그래프 아래 영역을 색칠하기
hold on;
rgb_light_green = [0.56, 0.93, 0.56]; % 연한 파란색 (R, G, B 값 순서)

area(x_seconds_shaded, shaded_values, 'FaceColor', rgb_light_green, 'FaceAlpha', 0.7);

%% 통계
n=nnz(systemWideAccessStatus)
systemWideAccessDuration=n*sc.SampleTime % seconds
scenarioDuration=seconds(sc.StopTime - sc.StartTime)
systemWideAccessPercentage=(systemWideAccessDuration/scenarioDuration)*100

%% play(sc)

orbitalElements(sat(1));


% 연속된 0의 가장 긴 길이 초기화
max_zero_length = 0;

% 현재 연속된 0의 길이 초기화
current_zero_length = 0;

% 행렬을 반복하면서 연속된 0의 가장 긴 길이를 찾음
for i = 1:length(systemWideAccessStatus)
    if systemWideAccessStatus(i) == 0
        current_zero_length = current_zero_length + 1;
        % 현재 연속된 0의 길이가 최대값보다 크다면 업데이트
        if current_zero_length > max_zero_length
            max_zero_length = current_zero_length;
        end
    else
        % 0이 아닌 값이 나올 때 현재 연속된 0의 길이 초기화
        current_zero_length = 0;
    end
end

max_zero_length*2


% 240열부터 480열까지의 부분 행렬 추출
subset = systemWideAccessStatus(240:480);

% 연속되는 1의 길이를 찾기 위한 변수 초기화
max_length = 0;
current_length = 0;

% 부분 행렬에서 연속되는 1의 최대 길이를 찾기 위한 반복문
for i = 1:length(subset)
    if subset(i) == 0
        current_length = current_length + 1;
        if current_length > max_length
            max_length = current_length;
        end
    else
        current_length = 0;
    end
end

max_length*2
