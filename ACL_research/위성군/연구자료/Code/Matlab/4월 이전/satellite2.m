startTime = datetime('2024-01-01 00:00:00');
stopTime = datetime('2024-01-02 00:00:00');
sampleTime = 120; % seconds
sc = satelliteScenario(startTime,stopTime,sampleTime)

%% 위성 정보

NP = 12
N = 22

for k = 0 : N-1

semiMajorAxis = 8059010;                                                                    % meters
eccentricity = 0;
inclination = 102.9;                                                                            % degrees
rightAscensionOfAscendingNode = round(720/ N *k,0)*0.5;                                                          % degrees
argumentOfPeriapsis = 98.3;                                                                    % degrees
trueAnomaly = 98.3 - 360 - NP * (360/N * 21);                                                                           % degrees
sat(k+1) = satellite(sc,semiMajorAxis,eccentricity,inclination,rightAscensionOfAscendingNode, ...
    argumentOfPeriapsis,trueAnomaly)

end



%%
% tleFile = "leoSatelliteConstellation.tle";
% sat = satellite(sc,tleFile)

names = sat.Name + " Camera";
cam = conicalSensor(sat,"Name",names,"MaxViewAngle",179.9)

name = "Ground Station";
latitude = 34.75; % 위도
longitude = 275.61; % 경도
altitude = 0; % 고도
minElevationAngle = 5; % degrees
geoSite = groundStation(sc, ...
    "Name",name,"Latitude", latitude, "Longitude", longitude, "Altitude", altitude, ...
     "MinElevationAngle",minElevationAngle);


ac = access(cam,geoSite);

% Properties of access analysis objects
ac(1);

v = satelliteScenarioViewer(sc,"ShowDetails",true);
sat(1).ShowLabel = true;
geoSite.ShowLabel = true;
show(sat(1));


fov = fieldOfView(cam([cam.Name] == "Satellite 1 Camera"))
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

% %% 타겟지향
% pointAt(sat,geoSite);
% 
% % Calculate system-wide access status
% for idx = 1:numel(ac)
%     [s,time] = accessStatus(ac(idx));
%     
%     if idx == 1
%         % Initialize system-wide access status vector in the first iteration
%         systemWideAccessStatus = s;
%     else
%         % Update system-wide access status vector by performing a logical OR
%         % with access status for the current camera-site combination
%         systemWideAccessStatus = or(systemWideAccessStatus,s);
%     end
% end
% 
% % Calculate system-wide access percentage
% n = nnz(systemWideAccessStatus);
% systemWideAccessDuration = n*sc.SampleTime;
% systemWideAccessPercentageWithTracking = (systemWideAccessDuration/scenarioDuration)*100

%% play(sc)


orbitalElements(sat(1));