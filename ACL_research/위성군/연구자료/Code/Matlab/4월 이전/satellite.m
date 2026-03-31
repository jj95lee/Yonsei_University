startTime = datetime('2024-01-01 00:00:00');
stopTime = datetime('2024-01-02 00:00:00');
sampleTime = 120; % seconds
sc = satelliteScenario(startTime,stopTime,sampleTime)

%% 위성 정보

semiMajorAxis = 7578137;                                                                    % meters
eccentricity = 0;
inclination = 102.9;                                                                           % degrees
rightAscensionOfAscendingNode = 0;                                                          % degrees
argumentOfPeriapsis = 98.3;                                                                    % degrees
trueAnomaly = 0;                                                                           % degrees
sat(1) = satellite(sc,semiMajorAxis,eccentricity,inclination,rightAscensionOfAscendingNode, ...
    argumentOfPeriapsis,trueAnomaly)

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

timeInSeconds = 1:721;  % 721개의 데이터 포인트, 86400초를 120초 간격으로 샘플링

% 원래의 timedate 값을 초 단위로 변환
time = seconds((timeInSeconds - 1) * sampleTime);  % 0을 포함한 720개의 값

plot(timeInSeconds, systemWideAccessStatus, "LineWidth", 2);
grid on;
xlabel("Time (seconds)");  % x축을 초 단위로 표시
ylabel("System-Wide Access Status");
hold on;

% 시간 단위로 변환하여 x축 눈금 설정
xticks(0:120:720);  % 120초 간격으로 x축 눈금 설정
xticklabels(string(0:120:720));  % x축 눈금 레이블을 초 단위로 표시


for i = 1 : 22
    plot(mod(timeInSeconds + round(720/22*i), 720), systemWideAccessStatus, "LineWidth", 2);
    grid on;
    xlabel("Time (seconds)");  % x축을 초 단위로 표시
    ylabel("System-Wide Access Status");
    hold on;
end

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