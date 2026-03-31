%% EpochTime

startTime = datetime('2024-01-01 00:00:00');
stopTime = datetime('2024-01-02 00:00:00');
sampleTime = 120; % seconds
sc = satelliteScenario(startTime,stopTime,sampleTime)

%% Orbital elements

semiMajorAxis = 8040137; % meters
eccentricity = 0;
inclination = 102.9; % degrees
rightAscensionOfAscendingNode = 0; % degrees
argumentOfPeriapsis = 98.3; % degrees
trueAnomaly = 0; % degrees
sat(1) = satellite(sc,semiMajorAxis,eccentricity,inclination,rightAscensionOfAscendingNode, ...
    argumentOfPeriapsis,trueAnomaly)

% tleFile = "leoSatelliteConstellation.tle";
% sat = satellite(sc,tleFile)

%% Sensor
names = sat.Name + " Camera";
cam = conicalSensor(sat,"Name",names,"MaxViewAngle",179.9)

%% Ground Station
name = "Ground Station";
latitude = 34.75; % latitude
longitude = 275.61; % longitude
altitude = 0; % altitude
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
% plot(x_seconds/120, systemWideAccessStatus, "LineWidth", 2);
% xlim([0 720]);
% 
% grid on;
% xlabel("Time (seconds)");  % 초 단위로 x축 표시
% ylabel("System-Wide Access Status");
% hold on;


% 아래 그래프를 색칠할 때 사용할 벡터 생성
x_seconds_shaded = x_seconds/120; % x축 벡터

% System-Wide Access Status 벡터를 0 또는 1 값으로 변환
shaded_values = systemWideAccessStatus; % 여기에 색칠할 영역의 값이 있어야 함

% 1인 값(접근 가능한 상태)을 0.5로 변환하여 색칠
shaded_values(shaded_values == 1) = 0.5;



d1 = [1	,2	,1	,1	,1	,1	,1	,1	,2	,2	,2	,1	,1	,2	,2	,2	,2	,2	,1	,2	,3	,2	,1	,1	,1	,1	,1	,2	,2	,2	,1	,1	,1	,2	,2	,1	,1	,1	,1	,1	,1	,2	,2	,2	,1	,1	,2	,2	,2	,2	,1	,1	,2	,2	,1	,1	,1	,1	,1	,2	,2	,2	,1	,1	,1	,1	,2	,1	,1	,1	,1	,1	,1	,2	,2	,2	,1	,1	,2	,2	,2	,2	,2	,1	,2	,2	,2	,1	,1	,1	,1	,1	,2	,2	,2	,1	,1	,1	,1	,2	,1	,1	,1	,1	,1	,1	,2	,2	,2	,1	,1	,2	,2	,2	,2	,2	,1	,2	,3	,2	,1	,1	,1	,1	,2	,2	,2	,2	,1	,1	,1	,2	,2	,1	,1	,1	,1	,1	,1	,2	,2	,2	,1	,1	,2	,2	,2	,2	,1	,2	,2	,2	,1	,1	,1	,1	,1	,2	,2	,2	,1	,1	,1	,1	,2	,1	,1	,1	,1	,1	,1	,2	,2	,2	,1	,1	,2	,2	,2	,2	,2	,1	,2	,3	,2	,1	,1	,1	,1	,1	,2	,2	,2	,1	,1	,1	,1	,2	,1	,1	,1	,1	,1	,1	,2	,2	,2	,1	,1	,2	,2	,2	,2	,1	,1	,2	,3	,2	,1	,1	,1	,1	,2	,2	,2	,1	,1	,1	,1	,2	,1	,1	,1	,1	,1	,1	,2	,2	,2	,2	,1	,2	,2	,2	,2	,2	,1	,2	,2	,2	,1	,1	,1	,1	,1	,2	,2	,2	,1	,1	,1	,1	,2	,1	,1	,1	,1	,1	,1	,2	,2	,2	,1	,1	,2	,2	,2	,2	,2	,1	,2	,3	,2	,1	,1	,1	,1	,1	,2	,2	,2	,1	,1	,1	,2	,2	,1	,1	,1	,1	,1	,1	,2	,2	,2	,1	,1	,2	,2	,2	,2	,1	,1	,2	,2	,2	,1	,1	,1	,1	,2	,2	,2	,1	,1	,1	,1	,2	,1	,1	,1	,1	,1	,1	,2	,2	,2	,1	,1	,2	,2	,2	,2	,2	,1	,2	,3	,2	,1	,1	,1	,1	,1	,2	,2	,2	,1	,1	,1	,1	,2	,1	,1	,1	,1	,1	,1	,2	,2	,2	,1	,1	,2	,2	,2	,2	,2	,1	,2	,3	,2	,1	,1	,1	,1	,2	,2	,2	,2	,1	,1	,1	,2	,2	,1	,1	,1	,1	,1	,1	,2	,2	,2	,1	,1	,2	,2	,2	,2	,1	,2	,2	,2	,1	,1	,1	,1	,1	,2	,2	,2	,1	,1	,1	,1	,2	,1	,1	,1	,1	,1	,1	,2	,2	,2	,1	,1	,2	,2	,2	,2	,2	,1	,2	,3	,2	,1	,1	,1	,1	,1	,2	,2	,2	,1	,1	,1	,1	,2	,1	,1	,1	,1	,1	,1	,2	,2	,2	,1	,1	,2	,2	,2	,2	,1	,1	,2	,3	,2	,1	,1	,1	,1	,2	,2	,2	,1	,1	,1	,1	,2	,1	,1	,1	,1	,1	,1	,2	,2	,2	,2	,1	,2	,2	,2	,2	,2	,1	,2	,2	,2	,1	,1	,1	,1	,1	,2	,2	,2	,1	,1	,1	,1	,2	,1	,1	,1	,1	,1	,1	,2	,2	,2	,1	,1	,2	,2	,2	,2	,2	,1	,2	,3	,2	,1	,1	,1	,1	,1	,2	,2	,2	,1	,1	,1	,2	,2	,1	,1	,1	,1	,1	,1	,2	,2	,2	,1	,1	,2	,2	,2	,2	,1	,1	,2	,2	,1	,1	,1	,1	,1	,2	,2	,2	,1	,1	,1	,1	,2	,1	,1	,1	,1	,1	,1	,2	,2	,2	,1	,1	,2	,2	,2	,2	,2	,1	,2	,2	,2	,1	,1	,1	,1	,1	,2	,2	,2	,1	,1	,1	,1	,2	,1	,1	,1	,1	,1	,1	,2	,2	,2	,1	,1	,2	,2	,2	,2	,2	,1	,2	,3	,2	,1	,1	,1	,1	,2	,2	,2	,2	,1	,1	,1	,2	,2	,1	,1	,1	,1	,1	,1	,2	,2	,2	,1	,1	,2	,2	,2	,2	,1	,2	,2	,2	,1	,1	,1	,1	,1	,2	,2	,2	,1	,1	,1	,1	,2	,1	,1	,1	,1	,1	,1	,2	,2	,2	,1	,1	,2	,2	,2	,2	,2	,1	,2	,2	,2	,1	,1	,1	,1	,1	,2	,2	,2	,1	,1	,1	,2	,2];
% a(22) = 687;
% b = systemWideAccessStatus(end-a(22)+1:end);
% c = systemWideAccessStatus(1:end-a(22));
% d22 = [b, c];
% 
% plot(x_seconds/120, d22, "LineWidth", 2);
% xlim([0 720]);




% 두 벡터를 더하여 합칩니다.


% 그래프를 그리기
plot(x_seconds/120, d1, "LineWidth", 3, "Color", '#77AC30');
xlim([0 720]);
ylim([0 3]);
yline(1, "--", "LineWidth", 2);
hold on;

grid on;
xlabel("n");  % 초 단위로 x축 표시
ylabel("bj[n]");
title('Coverage timeline');

% 그래프 아래 영역 색칠 및 투명도 조절
area(x_seconds/120, d1, 'FaceColor', '#E4EFD6', 'FaceAlpha', 0.5);


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
% Calculate system-wide access percentage

n = nnz(systemWideAccessStatus);
systemWideAccessDuration = n*sc.SampleTime;
systemWideAccessPercentageWithTracking = (systemWideAccessDuration/scenarioDuration)*100

%% play(sc)

orbitalElements(sat(1));

(720-14)/720*100
