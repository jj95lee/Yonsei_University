close all
clear
clc
str_filename1 = 'ReportFile_1day_LTAN22_20250101_SAT_SUN_ECI_221124.txt';
str_filename2 = 'ReportFile_1day_LTAN22_20250101_SAT_GS_ECEF_221124.txt';

%% Parsing GMAT data
delimiter = ';';
dataStartLine = 1;
opts = delimitedTextImportOptions('Delimiter', delimiter, 'DataLines', dataStartLine);
raw_scan1 = readmatrix(str_filename1,opts);
raw_scan2 = readmatrix(str_filename2,opts);
len = size(raw_scan1,1);
% 
% for ind = 1:len
% 	dataGMAT.UTCGregorian(ind,:) = raw_scan1{ind,1};
% 	dataGMAT.MJD_GMAT(ind,1)	 = str2double(raw_scan1{ind,2});
% 	dataGMAT.RV_ECI_X(ind,1)	 = str2double(raw_scan1{ind,3})*1e3;	% [km] -> [m]
% 	dataGMAT.RV_ECI_Y(ind,1)	 = str2double(raw_scan1{ind,4})*1e3;	% [km] -> [m]
% 	dataGMAT.RV_ECI_Z(ind,1)	 = str2double(raw_scan1{ind,5})*1e3;	% [km] -> [m]
% 	dataGMAT.RV_ECI_VX(ind,1)	 = str2double(raw_scan1{ind,6})*1e3;	% [km/s] -> [m/s]
% 	dataGMAT.RV_ECI_VY(ind,1)	 = str2double(raw_scan1{ind,7})*1e3;	% [km/s] -> [m/s]
% 	dataGMAT.RV_ECI_VZ(ind,1)	 = str2double(raw_scan1{ind,8})*1e3;	% [km/s] -> [m/s]
% 	dataGMAT.SUN_ECI_X(ind,1)	 = str2double(raw_scan1{ind,9})*1e3;	% [km] -> [m]
% 	dataGMAT.SUN_ECI_Y(ind,1)	 = str2double(raw_scan1{ind,10})*1e3;   % [km] -> [m]
% 	dataGMAT.SUN_ECI_Z(ind,1)	 = str2double(raw_scan1{ind,11})*1e3;   % [km] -> [m]
% 	dataGMAT.RV_ECEF_X(ind,1)	 = str2double(raw_scan2{ind,3})*1e3;	% [km] -> [m]
% 	dataGMAT.RV_ECEF_Y(ind,1)	 = str2double(raw_scan2{ind,4})*1e3;	% [km] -> [m]
% 	dataGMAT.RV_ECEF_Z(ind,1)	 = str2double(raw_scan2{ind,5})*1e3;	% [km] -> [m]
% 	dataGMAT.RV_ECEF_VX(ind,1)	 = str2double(raw_scan2{ind,6})*1e3;	% [km/s] -> [m/s]
% 	dataGMAT.RV_ECEF_VY(ind,1)	 = str2double(raw_scan2{ind,7})*1e3;	% [km/s] -> [m/s]
% 	dataGMAT.RV_ECEF_VZ(ind,1)	 = str2double(raw_scan2{ind,8})*1e3;	% [km/s] -> [m/s]
% 	dataGMAT.GSYS_ECEF_X(ind,1)	 = str2double(raw_scan2{ind,9})*1e3;	% [km] -> [m]
% 	dataGMAT.GSYS_ECEF_Y(ind,1)	 = str2double(raw_scan2{ind,10})*1e3;   % [km] -> [m]
% 	dataGMAT.GSYS_ECEF_Z(ind,1)	 = str2double(raw_scan2{ind,11})*1e3;   % [km] -> [m]
% end
% 
% dataGMAT.ElapsedDays = dataGMAT.MJD_GMAT - dataGMAT.MJD_GMAT(1);
% dataGMAT.ElapsedSecs = dataGMAT.ElapsedDays*86400;






% %% Vector in ECI
% pos_SAT_ECI = [dataGMAT.RV_ECI_X,	dataGMAT.RV_ECI_Y,	dataGMAT.RV_ECI_Z];
% pos_SUN_ECI = [dataGMAT.SUN_ECI_X,	dataGMAT.SUN_ECI_Y, dataGMAT.SUN_ECI_Z];
% pos_rel_SAT2SUN_ECI = pos_SUN_ECI - pos_SAT_ECI;
% 
% %% Vector in ECEF
% pos_SAT_ECEF = [dataGMAT.RV_ECEF_X,	dataGMAT.RV_ECEF_Y,	dataGMAT.RV_ECEF_Z];
% pos_GSYS_ECEF = [dataGMAT.GSYS_ECEF_X,	dataGMAT.GSYS_ECEF_Y, dataGMAT.GSYS_ECEF_Z];
% pos_rel_GSYS2SAT_ECEF = pos_SAT_ECEF - pos_GSYS_ECEF;
% range_GSYS2SAT = sqrt(sum(pos_rel_GSYS2SAT_ECEF.^2,2));
% 
% %% Get Eclipse
% tmp1_n = sqrt(sum(pos_rel_SAT2SUN_ECI.^2,2));
% tmp2_n = sqrt(sum(pos_SAT_ECI.^2,2));
% ang_SAT2SUN = acos(sum(pos_rel_SAT2SUN_ECI.*pos_SAT_ECI,2)./tmp1_n./tmp2_n);
% tmp_RE	= 6378000; % [m] tmp_RE > tmp_RE_GSYS & tmp_RE_GSHS; + uncertainty
% tmp_SAT = 6878000; % [m] SMA SAT
% alpha = acos(tmp_RE/tmp_SAT);
% flag_eclipse_SAT = ang_SAT2SUN > pi/2 + alpha;
% flag_dayside_SAT = ~flag_eclipse_SAT;
% 
% 
% %% index 8680 -> Time 86400
% 
%  x = 1;
%   for j = 1 : 8679
%      for i = 0 : dataGMAT.ElapsedSecs(j+1) - dataGMAT.ElapsedSecs(j) - 1
%         dataGMAT.data(x,1) = ang_SAT2SUN(j) + ((ang_SAT2SUN(j+1) - ang_SAT2SUN(j)) / (dataGMAT.ElapsedSecs(j+1) - dataGMAT.ElapsedSecs(j)))*i;
%         x = x+1;
%      end
%   end
% 
% %% Make Topocentric Coordinate System (Ehat == East direction / Nhat == North direction / nhat == Zenith direction)
% % @GSYS
% RE_a = 6378137;
% RE_b = 6356752.314245;
% pos_GSYS_SPH = [ 37.563958, 126.93475, 90];
% pos_rel_GSYS2SAT_RAE = zeros(len,3);
% R_ECEF2GSYS = [0, 1, 0;0, 0, 1;1, 0, 0]*Rot2(-pos_GSYS_SPH(1)*pi/180)*Rot3(pos_GSYS_SPH(2)*pi/180);
% pos_rel_GSYS2SAT_ENn = pos_rel_GSYS2SAT_ECEF*R_ECEF2GSYS';
% 
% for ind = 1:len
% 	pos_GSYS_ECEF_k = pos_GSYS_ECEF(ind,:);
% 	pos_rel_GSYS2SAT_ECEF_k = pos_rel_GSYS2SAT_ECEF(ind,:);
% 	pos_rel_GSYS2SAT_ENn_k = (R_ECEF2GSYS * pos_rel_GSYS2SAT_ECEF_k')';
% 	pos_rel_GSYS2SAT_ENn(ind,:) = pos_rel_GSYS2SAT_ENn_k;
% 	rho_rel_GSYS2SAT_k = norm(pos_rel_GSYS2SAT_ENn_k);
% 	elevation_rel_GSYS2SAT_k = asind(pos_rel_GSYS2SAT_ENn_k(3)/rho_rel_GSYS2SAT_k);
% 	gndposSAT_tmp1 = [pos_rel_GSYS2SAT_ENn_k(1,1:2), 0];
% 	azimuth_rel_GSYS2SAT_k = acosd(dot([0,1,0],gndposSAT_tmp1)/norm(gndposSAT_tmp1));
% 	if gndposSAT_tmp1(1) < 0
% 		azimuth_rel_GSYS2SAT_k = 360 - azimuth_rel_GSYS2SAT_k;
% 	end
% 	pos_rel_GSYS2SAT_RAE_k = [rho_rel_GSYS2SAT_k, azimuth_rel_GSYS2SAT_k, elevation_rel_GSYS2SAT_k];
% 	pos_rel_GSYS2SAT_RAE(ind,:) = pos_rel_GSYS2SAT_RAE_k;
% end
% 
% %% Divide Mission Period
% % Find contact pass
% flag_el_pos = pos_rel_GSYS2SAT_RAE(:,3) > 0;
% ind_str_el_pos = find(diff(flag_el_pos) == 1);
% ind_end_el_pos = find(diff(flag_el_pos) == -1);
% k = 1;
% 
% for ind = 1:length(ind_str_el_pos)
% 	if max(pos_rel_GSYS2SAT_RAE(ind_str_el_pos(ind):ind_end_el_pos(ind),3)) > 30
% 		ind_str_contact(k,1) = ind_str_el_pos(ind);
% 		ind_end_contact(k,1) = ind_end_el_pos(ind);
% 		k = k + 1;
% 	end
% end
% 
% flag_contact = false(len,1);
% for ind = 1:k-1
% 	flag_contact(ind_str_contact(ind):ind_end_contact(ind),1) = true;
% end
% 
% % Sun pointing mode
% flag_sunpointing = ~flag_contact & flag_dayside_SAT;
% 
% % Nadir pointing mode
% flag_nadirpointing = ~flag_sunpointing;
% flag_mode = flag_sunpointing + flag_nadirpointing * 2; % 1: sun pointing, 2: nadir pointing
% 
% %% Get Angle b.t.w. Sun & Solar Cells
% % 0		rad -> full charging / pointing Sun
% % pi/2	rad -> non-charging
% % @Eclipse:: all values are "pi/2"
% %--Basic Equation for Nadir Pointing--%
% % angle_SUN_CELL = acos(dot(pos_SAT_ECI,pos_rel_SAT2SUN_ECI)/norm(pos_SAT_ECI)/norm(pos_rel_SAT2SUN_ECI));
% %-------------------------------------%
% angle_SUN_CELL = acos(sum(pos_SAT_ECI.*pos_rel_SAT2SUN_ECI,2) ./ sqrt(sum(pos_SAT_ECI.^2,2)) ./ sqrt(sum(pos_rel_SAT2SUN_ECI.^2,2))); % [rad]
% 
% %angle_SUN_CELL(angle_SUN_CELL > pi/2,1) = pi/2;	% Cut values over pi/2
% %angle_SUN_CELL(flag_eclipse_SAT,1) = pi/2;		% Set to pi/2 @Eclipse
% %angle_SUN_CELL(flag_sunpointing,1) = 0;			% Set to 0 @Sun pointing
% 
% %% Plot Results
% f0 = figure();
% set(f0,'defaultAxesColorOrder',[0,0,0;0,0,1]);
% 
% yyaxis left
% plot(dataGMAT.ElapsedSecs, angle_SUN_CELL, 'LineWidth', 1.2)
% grid on
% xlabel('time [sec]')
% ylabel('Angle b.t.w. Sun & Cells')
% title('Data for 24 hours / LTAN22h')
% 
% yyaxis right
% plot(dataGMAT.ElapsedSecs,	flag_mode, 'LineStyle','--')
% ylabel('mode')
% ylim([0 2.5])
% yticks([1 2])
% 
% %% Plot Results
% f1 = figure();
% set(f1,'defaultAxesColorOrder',[0,0,0;0,0,1]);
% 
% yyaxis left
% plot(angle_SUN_CELL, 'LineWidth', 1.2)
% grid on
% xlabel('time [sec]')
% ylabel('Angle b.t.w. Sun & Cells')
% title('Data for 24 hours / LTAN22h')
% 
% yyaxis right
% plot(flag_mode, 'LineStyle','--')
% ylabel('mode')
% ylim([0 2.5])
% yticks([1 2])