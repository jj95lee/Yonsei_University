clear all
clc

NP = 12
N = 22


storedValues = zeros(1, 22);


for k = 1:22
                                                                   % degrees
trueAnomaly = 98.3 + 360 - NP * (360/N * k);  

storedValues(k) = trueAnomaly;

end

% 저장된 값 출력
disp('저장된 값들:');
disp(storedValues);

