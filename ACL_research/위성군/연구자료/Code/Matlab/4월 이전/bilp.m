clear all
clc

startTime = datetime('2024-01-01 00:00:00');
stopTime = datetime('2024-01-02 00:00:00');
sampleTime = 120; % seconds
sc = satelliteScenario(startTime,stopTime,sampleTime)

%% 위성 정보

NP = 14;
N = 29;

L = 720;

for k = 0:N-1
semiMajorAxis = 7274140;                                                                    % meters
eccentricity = 0;
inclination = 50;                                                                            % degrees                                                       % degrees
rightAscensionOfAscendingNode = round(L/N)*k; 
argumentOfPeriapsis = 0;                                                                    % degrees                                                                         % degrees
trueAnomaly = mod(-NP*rightAscensionOfAscendingNode,720);
sat(k+1) = satellite(sc,semiMajorAxis,eccentricity,inclination,rightAscensionOfAscendingNode, ...
    argumentOfPeriapsis,trueAnomaly,"OrbitPropagator","two-body-keplerian")

 % 작업공간에 저장
 RAAN(k+1) = rightAscensionOfAscendingNode;   
 Ta(k+1) = trueAnomaly;
 n_pos(k+1) = round(L/N*k,0); 

end

% semiMajorAxis = 8048280;                                                                    % meters
% eccentricity = 0;
% inclination = 50;                                                                           % degrees
% rightAscensionOfAscendingNode = 0;                                                          % degrees
% argumentOfPeriapsis = 0;                                                                    % degrees
% trueAnomaly = 96;                                                                           % degrees
% sat(1) = satellite(sc,semiMajorAxis,eccentricity,inclination,rightAscensionOfAscendingNode, ...
%     argumentOfPeriapsis,trueAnomaly)
% 
% semiMajorAxis = 8048280;                                                                    % meters
% eccentricity = 0;
% inclination = 50;                                                                           % degrees
% rightAscensionOfAscendingNode = 33*0.5;                                                          % degrees
% argumentOfPeriapsis = 0;                                                                    % degrees
% trueAnomaly = 264;                                                                           % degrees
% sat(2) = satellite(sc,semiMajorAxis,eccentricity,inclination,rightAscensionOfAscendingNode, ...
%     argumentOfPeriapsis,trueAnomaly)
% 
% semiMajorAxis = 8048280;                                                                    % meters
% eccentricity = 0;
% inclination = 50;                                                                           % degrees
% rightAscensionOfAscendingNode = 65*0.5;                                                          % degrees
% argumentOfPeriapsis = 0;                                                                    % degrees
% trueAnomaly = 66;                                                                           % degrees
% sat(3) = satellite(sc,semiMajorAxis,eccentricity,inclination,rightAscensionOfAscendingNode, ...
%     argumentOfPeriapsis,trueAnomaly)
% 
% semiMajorAxis = 8048280;                                                                    % meters
% eccentricity = 0;
% inclination = 50;                                                                           % degrees
% rightAscensionOfAscendingNode = 98*0.5;                                                          % degrees
% argumentOfPeriapsis = 0;                                                                    % degrees
% trueAnomaly = 228;                                                                           % degrees
% sat(4) = satellite(sc,semiMajorAxis,eccentricity,inclination,rightAscensionOfAscendingNode, ...
%     argumentOfPeriapsis,trueAnomaly)
% 
% semiMajorAxis = 8048280;                                                                    % meters
% eccentricity = 0;
% inclination = 50;                                                                           % degrees
% rightAscensionOfAscendingNode = 131*0.5;                                                          % degrees
% argumentOfPeriapsis = 0;                                                                    % degrees
% trueAnomaly = 30;                                                                           % degrees
% sat(5) = satellite(sc,semiMajorAxis,eccentricity,inclination,rightAscensionOfAscendingNode, ...
%     argumentOfPeriapsis,trueAnomaly)
% 
% semiMajorAxis = 8048280;                                                                    % meters
% eccentricity = 0;
% inclination = 50;                                                                           % degrees
% rightAscensionOfAscendingNode = 164*0.5;                                                          % degrees
% argumentOfPeriapsis = 0;                                                                    % degrees
% trueAnomaly = 198;                                                                           % degrees
% sat(6) = satellite(sc,semiMajorAxis,eccentricity,inclination,rightAscensionOfAscendingNode, ...
%     argumentOfPeriapsis,trueAnomaly)
% 
% semiMajorAxis = 8048280;                                                                    % meters
% eccentricity = 0;
% inclination = 50;                                                                           % degrees
% rightAscensionOfAscendingNode = 196*0.5;                                                          % degrees
% argumentOfPeriapsis = 0;                                                                    % degrees
% trueAnomaly = 0;                                                                           % degrees
% sat(7) = satellite(sc,semiMajorAxis,eccentricity,inclination,rightAscensionOfAscendingNode, ...
%     argumentOfPeriapsis,trueAnomaly)
% 
% semiMajorAxis = 8048280;                                                                    % meters
% eccentricity = 0;
% inclination = 50;                                                                           % degrees
% rightAscensionOfAscendingNode = 229*0.5;                                                          % degrees
% argumentOfPeriapsis = 0;                                                                    % degrees
% trueAnomaly = 162;                                                                           % degrees
% sat(8) = satellite(sc,semiMajorAxis,eccentricity,inclination,rightAscensionOfAscendingNode, ...
%     argumentOfPeriapsis,trueAnomaly)
% 
% semiMajorAxis = 8048280;                                                                    % meters
% eccentricity = 0;
% inclination = 50;                                                                           % degrees
% rightAscensionOfAscendingNode = 262*0.5;                                                          % degrees
% argumentOfPeriapsis = 0;                                                                    % degrees
% trueAnomaly = 330;                                                                           % degrees
% sat(9) = satellite(sc,semiMajorAxis,eccentricity,inclination,rightAscensionOfAscendingNode, ...
%     argumentOfPeriapsis,trueAnomaly)
% 
% semiMajorAxis = 8048280;                                                                    % meters
% eccentricity = 0;
% inclination = 50;                                                                           % degrees
% rightAscensionOfAscendingNode = 295*0.5;                                                          % degrees
% argumentOfPeriapsis = 0;                                                                    % degrees
% trueAnomaly = 132;                                                                           % degrees
% sat(10) = satellite(sc,semiMajorAxis,eccentricity,inclination,rightAscensionOfAscendingNode, ...
%     argumentOfPeriapsis,trueAnomaly)
% 
% semiMajorAxis = 8048280;                                                                    % meters
% eccentricity = 0;
% inclination = 50;                                                                           % degrees
% rightAscensionOfAscendingNode = 327*0.5;                                                          % degrees
% argumentOfPeriapsis = 0;                                                                    % degrees
% trueAnomaly = 294;                                                                           % degrees
% sat(11) = satellite(sc,semiMajorAxis,eccentricity,inclination,rightAscensionOfAscendingNode, ...
%     argumentOfPeriapsis,trueAnomaly)
% 
% semiMajorAxis = 8048280;                                                                    % meters
% eccentricity = 0;
% inclination = 50;                                                                           % degrees
% rightAscensionOfAscendingNode = 360*0.5;                                                          % degrees
% argumentOfPeriapsis = 0;                                                                    % degrees
% trueAnomaly = 96;                                                                           % degrees
% sat(12) = satellite(sc,semiMajorAxis,eccentricity,inclination,rightAscensionOfAscendingNode, ...
%     argumentOfPeriapsis,trueAnomaly)
% 
% semiMajorAxis = 8048280;                                                                    % meters
% eccentricity = 0;
% inclination = 50;                                                                           % degrees
% rightAscensionOfAscendingNode = 393*0.5;                                                          % degrees
% argumentOfPeriapsis = 0;                                                                    % degrees
% trueAnomaly = 264;                                                                           % degrees
% sat(13) = satellite(sc,semiMajorAxis,eccentricity,inclination,rightAscensionOfAscendingNode, ...
%     argumentOfPeriapsis,trueAnomaly)
% 
% semiMajorAxis = 8048280;                                                                    % meters
% eccentricity = 0;
% inclination = 50;                                                                           % degrees
% rightAscensionOfAscendingNode = 425*0.5;                                                          % degrees
% argumentOfPeriapsis = 0;                                                                    % degrees
% trueAnomaly = 66;                                                                           % degrees
% sat(14) = satellite(sc,semiMajorAxis,eccentricity,inclination,rightAscensionOfAscendingNode, ...
%     argumentOfPeriapsis,trueAnomaly)
% 
% semiMajorAxis = 8048280;                                                                    % meters
% eccentricity = 0;
% inclination = 50;                                                                           % degrees
% rightAscensionOfAscendingNode = 458*0.5;                                                          % degrees
% argumentOfPeriapsis = 0;                                                                    % degrees
% trueAnomaly = 228;                                                                           % degrees
% sat(15) = satellite(sc,semiMajorAxis,eccentricity,inclination,rightAscensionOfAscendingNode, ...
%     argumentOfPeriapsis,trueAnomaly)
% 
% semiMajorAxis = 8048280;                                                                    % meters
% eccentricity = 0;
% inclination = 50;                                                                           % degrees
% rightAscensionOfAscendingNode = 491*0.5;                                                          % degrees
% argumentOfPeriapsis = 0;                                                                    % degrees
% trueAnomaly = 30;                                                                           % degrees
% sat(16) = satellite(sc,semiMajorAxis,eccentricity,inclination,rightAscensionOfAscendingNode, ...
%     argumentOfPeriapsis,trueAnomaly)
% 
% semiMajorAxis = 8048280;                                                                    % meters
% eccentricity = 0;
% inclination = 50;                                                                           % degrees
% rightAscensionOfAscendingNode = 524*0.5;                                                          % degrees
% argumentOfPeriapsis = 0;                                                                    % degrees
% trueAnomaly = 198;                                                                           % degrees
% sat(17) = satellite(sc,semiMajorAxis,eccentricity,inclination,rightAscensionOfAscendingNode, ...
%     argumentOfPeriapsis,trueAnomaly)
% 
% semiMajorAxis = 8048280;                                                                    % meters
% eccentricity = 0;
% inclination = 50;                                                                           % degrees
% rightAscensionOfAscendingNode = 556*0.5;                                                          % degrees
% argumentOfPeriapsis = 0;                                                                    % degrees
% trueAnomaly = 0;                                                                           % degrees
% sat(18) = satellite(sc,semiMajorAxis,eccentricity,inclination,rightAscensionOfAscendingNode, ...
%     argumentOfPeriapsis,trueAnomaly)
% 
% semiMajorAxis = 8048280;                                                                    % meters
% eccentricity = 0;
% inclination = 50;                                                                           % degrees
% rightAscensionOfAscendingNode = 589*0.5;                                                          % degrees
% argumentOfPeriapsis = 0;                                                                    % degrees
% trueAnomaly = 162;                                                                           % degrees
% sat(19) = satellite(sc,semiMajorAxis,eccentricity,inclination,rightAscensionOfAscendingNode, ...
%     argumentOfPeriapsis,trueAnomaly)
% 
% semiMajorAxis = 8048280;                                                                    % meters
% eccentricity = 0;
% inclination = 50;                                                                           % degrees
% rightAscensionOfAscendingNode = 622*0.5;                                                          % degrees
% argumentOfPeriapsis = 0;                                                                    % degrees
% trueAnomaly = 330;                                                                           % degrees
% sat(20) = satellite(sc,semiMajorAxis,eccentricity,inclination,rightAscensionOfAscendingNode, ...
%     argumentOfPeriapsis,trueAnomaly)
% 
% semiMajorAxis = 8048280;                                                                    % meters
% eccentricity = 0;
% inclination = 50;                                                                           % degrees
% rightAscensionOfAscendingNode = 655*0.5;                                                          % degrees
% argumentOfPeriapsis = 0;                                                                    % degrees
% trueAnomaly = 132;                                                                           % degrees
% sat(21) = satellite(sc,semiMajorAxis,eccentricity,inclination,rightAscensionOfAscendingNode, ...
%     argumentOfPeriapsis,trueAnomaly)

% semiMajorAxis = 8048280;                                                                    % meters
% eccentricity = 0;
% inclination = 50;                                                                           % degrees
% rightAscensionOfAscendingNode = 687*0.5;                                                          % degrees
% argumentOfPeriapsis = 0;                                                                    % degrees
% trueAnomaly = 294;                                                                           % degrees
% sat(22) = satellite(sc,semiMajorAxis,eccentricity,inclination,rightAscensionOfAscendingNode, ...
%     argumentOfPeriapsis,trueAnomaly)




%%
% tleFile = "leoSatelliteConstellation.tle";
% sat = satellite(sc,tleFile)

names = sat.Name + " Camera";
cam = conicalSensor(sat,"Name",names,"MaxViewAngle",179.9999999999999)

name = "Ground Station";
latitude = 39.0392; % 위도
longitude = 125.7625; % 경도
altitude = 0; % 고도
minElevationAngle = 10; % degrees
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
figure;
plot(x_seconds/120-100, systemWideAccessStatus, "LineWidth", 2);
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

%% V

shiftAmount = 1;

V = [];

for i = 1:length(systemWideAccessStatus)
    newRow = circshift(systemWideAccessStatus, [0, i-shiftAmount]);
    V = [V; newRow];
end


%% BILP position

% 행렬 크기 정의
row_count = 720;
col_count = 1;

% 모두 0으로 채워진 행렬 생성
Z = zeros(row_count, col_count);

% 덮어쓸 특정 행 지정
% target_rows = [39, 73, 79, 89, 170, 184, 234, 250, 331, 341, 347, 492, 502, 542, 638, 648, 654, 663];
target_rows = [18, 114, 124, 130, 139, 235, 269, 275, 285, 366, 380, 430, 446, 527, 537, 543, 688, 698];

% 특정 행에 1 덮어쓰기
Z(target_rows, :) = 1;

%% Binary random
numOnes = 18;
numZeros = 720-numOnes;
count = 0;

% while true
% 
% W = [zeros(numZeros, 1); ones(numOnes, 1)];
% 
% W = W(randperm(length(W)));

X = V*Z; 
Y = sum(X >= 1)/720*100;

%     if Y >= 95
%         break;
%     end
% 
%     count = count + 1;
% end

disp(['Y% 이상이 나올 때까지의 반복 횟수: ' num2str(count)]);

x_values = 1:720;
y_values = X(1:720);

% 그래프 다시 그리기
figure;
plot(x_values, y_values, "LineWidth", 2);
xlim([0 720]);
ylim([0 4]);
grid on;
xlabel("n");  % 초 단위로 x축 표시
ylabel("Revisit Time");
title('Total Coverage Timeline')

% 주어진 A 행렬
% A = [1, 1, 1, 0, 0, 1, 1, 0, 0, 1;
%      1, 1, 1, 1, 0, 0, 1, 1, 0, 0;
%      0, 1, 1, 1, 1, 0, 0, 1, 1, 0;
%      0, 0, 1, 1, 1, 1, 0, 0, 1, 1;
%      1, 0, 0, 1, 1, 1, 1, 0, 0, 1;
%      1, 1, 0, 0, 1, 1, 1, 1, 0, 0;
%      0, 1, 1, 0, 0, 1, 1, 1, 1, 0;
%      0, 0, 1, 1, 0, 0, 1, 1, 1, 1;
%      1, 0, 0, 1, 1, 0, 0, 1, 1, 1;
%      1, 1, 0, 0, 1, 1, 0, 0, 1, 1];

% 원하는 1의 개수 설정
desired_ones = 18;

% 임계값 계산
threshold = desired_ones / 720;

% % C 행렬의 조건과 1의 개수를 충족하는 B 행렬 찾기
% C_condition = false;
% 
% while ~C_condition
%     % 0의 확률을 높게 설정
%     B = rand(720, 1) < threshold;
%     while sum(B) > desired_ones % 1의 개수가 N 이하일 때까지 반복
%         B = rand(720, 1) < threshold;
%     end
%     C = V * B; % C 행렬 계산
%     C_condition = all(sum(C(:) >= 1)/720 >= 10); % C 행렬의 조건 확인


% B 행렬 그래프
subplot(2, 1, 1);
plot(B, '-o r', 'LineWidth', 1, 'MarkerSize', 1);
title('Constellation Pattern Vector');
xlim([0, 720]);
ylim([-0.2, 1.5]);

% C 행렬 그래프
subplot(2, 1, 2);
plot(find(C), C(C~=0), '-o', 'Color', [0, 0.5, 0], 'LineWidth', 1, 'MarkerSize', 1);
title('Coverage Timeline');
xlim([0, 720]);
ylim([0, 10]);

xlabel('n');
ylabel('');
