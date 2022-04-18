% clear all;clc;
c = physconst('LightSpeed');        % propagation speed
fc = 2.4e9;      % carrier frequency
lambda = c/fc;  % wavelength
antenna = phased.IsotropicAntennaElement(...
    'FrequencyRange',[3e8 3e9]);
txarray = phased.ULA('NumElements',2,'ElementSpacing',lambda/2);    
% viewArray(txarray);
txmipos_realPos = getElementPosition(txarray);
txmipos = txmipos_realPos/lambda;  % x, y, z, the default coordinate that used is y
% the reason divide by lambda is normalized 

pattern(txarray,fc,-180:180,0,'PropagationSpeed',c,...
    'CoordinateSystem','rectangular',...
    'Type','powerdb','Normalize',true)

%% array multiplificaiton

%%
freq = 20e6;
azangles = -180:180;
response = txarray(freq,azangles);
pattern(antenna,fc,azangles,[-90:0.1:90],...
    'Type','efield',...
    'CoordinateSystem','polar')
% patternAzimuth(antenna,fc,30) % elevation = 0
% patternAzimuth(antenna,fc,30,'Azimuth',[-20:20]) % % elevation = 30
% patternElevation(antenna,fc,45)
%% signal delay between array element
delay = phased.ElementDelay('SensorArray',txarray);
tau = delay([-90;0]); % azimuth = -90, elevation=0     du: -> d*sin(theta)/c = -2.0834e-10

%%
% If IncludeElementResponse is false, the computation assumes that the elements are isotropic and SV does not include the individual element responses.
txarraystv = phased.SteeringVector('SensorArray',txarray,'PropagationSpeed',c); 
txang = [80;0];  % tx angle (azimuth; elevation)          -180<=azimuth<=180  -90<=elevation<=90
wt = txarraystv(fc,txang)'; % steering vector, N-by-M-by-L, N is the number of array elements, M is the number of steering directions, and L is the number of frequency.

% txbeam_ang = 0;
% % txsv is an N-by-M complex -valued matrix, N is the number of element positions(antennas), M is the number of incoming waves
% txsv = steervec(txmipos,txbeam_ang);  % returns the sv for incoming plane wave impinging on a sensor array
% txbeam = wt * txsv
%% test 
wt_scale1 = int16((wt * 2^7)-1);
wt_scale2 = (wt * 2^7);


%%
txbeam_ang = -90:90;
%txbeam_ang = -180:180;
txbeam_ang_rad = (pi*txbeam_ang)/180;
txsv = steervec(txmipos(2,:),txbeam_ang);
txbeam1 = abs(wt*txsv);   %% wt_scale wt
txbeam2 = abs(wt*txsv);   %% wt_scale wt

%txbeam = txbeam/max(txbeam);
%[txbeampos_x,txbeampos_y] = pol2cart(deg2rad(txbeam_ang),txbeam);

% hFig = figure(1);
% set(hFig, 'Position', [0 0 800 400]);
figure()
% subplot(1,2,1);
plot(txbeam_ang,txbeam1);
hold on
plot(txbeam_ang,txbeam2);
xlabel('Tx Beam Angle');ylabel('Normalized Amplitude');
set(gca,'xtick',-180:10:180)
%xlim([txbeam_ang(1) txbeam_ang(end)]); ylim([0 1.0]);
title('30 degree beam pattern')

wt = wt.'
% subplot(1,2,2);  
%    polarplot(txbeam_ang_rad,txbeam,'r');
%    set(gca,'RTickLabels',[]);










%% 绘制beam图
% x = -90:10:90;
% data_0_l = [78.563, 79.523, 80.187,  80.986, 82.023, 84.790, 85.541, 87.908, 88.565, 86.146]; %左0-90
% data_0_r = [80.097, 81.351, 82.572,  84.613, 85.207, 85.912, 86.937, 85.445, 84.932];
% data_0 = [flip(data_0_l) data_0_r];
% figure(2)
% subplot(1,2,1);
% plot(x,(data_0/(max(data_0)))*-1 +1.88);
% title('2 Antennas with 0 degree');
% xlabel('Tx beam angle');ylabel('Normalized power');
% % values = spcrv([[x(1) x x(end)];[data_0(1) data_0 data_0(end)]],2);
% % plot(values(1,:),values(2,:));
% 
% %%30度
% xx = -90:5:90;
% data_30_l = [84.091, 84.357, 85.133, 85.733, 86.576, 87.566, 89.553, 88.968, 88.365, 87.765, 86.751, 84.413, 83.269, 83.341, 84.131, 85.432, 86.235, 87.167, 87.266];
% data_30_r = [83.844, 83.504, 83.137, 82.451, 82.167, 81.779, 82.624, 82.866, 83.223, 83.518, 84.156, 84.468, 85.465, 86.090, 86.788, 87.564, 88.593, 89.167];
% data_30 =  [flip(data_30_l) data_30_r];
% subplot(1,2,2);
% plot(xx,(data_30/max(data_30))*-1 + 1.91);
% title('2 Antennas with 30 degree');
% xlabel('Tx beam angle');ylabel('Normalized power');
% %(data_30/(-1*min(data_30)))
% % figure(4)
% % values = spcrv([[xx(1) xx xx(end)];[data_30(1) data_30 data_30(end)]],2);
% % plot(values(1,:),values(2,:));

% %%4 antennas
% %%0 degree
% xxx = -90:5:90;
% data4_0_l = [70.567, 70.483, 71.031, 71.863, 72.446, 73.536, 74.613, 75.545, 76.247, 77.995, 78.467, 79.257, 80.334, 81.275, 81.868, 82.356, 82.886, 82.411, 81.877];
% data4_0_r = [70.686, 71.015, 71.667, 72.441, 73.013, 73.996, 74.795, 75.531, 76.232, 76.976, 77.406, 78.576, 79.484, 80.541, 81.011, 81.786, 82.231, 81.657];
% data4_0 =  [flip(data4_0_l) data4_0_r];
% figure(4)
% subplot(1,2,1);
% plot(xxx,(data4_0/max(data4_0))*-1 +1.85);
% %(data4_0/max(data4_0))*-1
% title('4 Antennas with 0 degree');
% xlabel('Tx beam angle');ylabel('Normalized power');
% %%30度
% xx = -90:5:90;
% data4_30_l = [67.831, 68.950, 69.434, 70.853, 71.431, 70.024, 69.899, 67.737, 66.772, 67.677, 68.551, 69.351, 69.971, 70.451, 71.257, 73.629, 74.445, 75.541, 77.519];
% data4_30_r = [66.721, 65.957, 64.731, 62.995, 61.434, 60.850, 61.275, 61.945, 62.741, 63.899, 64.422, 65.344, 66.477, 68.895, 69.970, 72.713, 74.843, 76.354];
% data4_30 =  [flip(data4_30_l) data4_30_r];
% subplot(1,2,2);
% plot(xx,(data4_30/max(data4_30))*-1 + 1.76);
% title('4 Antennas with 30 degree');
% xlabel('Tx beam angle');ylabel('Normalized power');

%%0 degree
% xxx = -90:10:90;
% data4_0_l = [72.168, 73.544, 75.269, 77.441, 78.541, 77.644, 77.241, 77.544, 76.441, 71.647];
% data4_0_r = [72.641, 73.441, 72.841, 71.644, 71.669, 69.699, 69.344, 68.477, 68.221];
% data4_0 =  [data4_0_l data4_0_r];
% subplot(1,2,2);
% plot(xxx,(data4_0/max(data4_0))*-1 +1.85);
% %(data4_0/max(data4_0))*-1
% title('4 Antennas with 0 degree');
% xlabel('Tx beam angle');ylabel('Normalized power');

% xxx = -90:10:90;
% data4_10_l = [76.144, 74.647, 72.641, 72.201, 71.241, 72.241, 70.345, 69.047, 66.261, 64.341];
% data4_10_r = [64.764, 65.041, 66.769, 67.362, 69.034, 69.544, 70.341, 70.697, 70.644];
% data4_10 =  [data4_10_l data4_10_r];
% figure(5)
% plot(xxx,(data4_10/max(data4_10))*-1 +1.83);
% title('4 Antennas with 0 degree');
% xlabel('Tx beam angle');ylabel('Normalized power');

% %40-50 degree
% xxx = -90:10:90;
% data4_10_l = [69.677, 71.369, 72.644, 74.761, 75.241, 75.041, 74.609, 73.441, 73.021, 71.691];
% data4_10_r = [70.745, 69.987, 67.241, 64.271, 63.477, 64.607, 64.871, 65.874, 65.969];
% data4_10 =  [data4_10_l data4_10_r];
% subplot(1,2,2);
% plot(xxx,(data4_10/max(data4_10))*-1 +1.83);
% title('4 Antennas with 50 degree');
% xlabel('Tx beam angle');ylabel('Normalized power');

%40 degree
% xxx = -90:10:90;
% data4_10_l = [68.055, 68.341, 68.861, 69.267, 68.571, 68.041, 67.799, 66.541, 65.977, 65.541];
% data4_10_r = [65.011, 64.744, 64.341, 64.027, 64.604, 66.266, 67.041, 67.631, 68.131];
% data4_10 =  [data4_10_l data4_10_r];
% subplot(1,2,2);
% plot(xxx,(data4_10/max(data4_10))*-1 +1.92);
% title('4 Antennas with 40 degree');
% xlabel('Tx beam angle');ylabel('Normalized power');

%%usrp 30度
% xxx = -90:10:90;
% data4_0_l = [58.649, 57.877, 57.017, 56.641, 56.322, 56.801, 57.901, 56.449, 57.471, 58.609];
% data4_0_r = [56.241, 55.407, 53.051, 54.407, 55.571, 56.441, 57.071, 57.631, 57.841];
% data4_0 =  [data4_0_l data4_0_r];
% subplot(1,2,2);
% plot(xxx,(data4_0/max(data4_0))*-1 +1.9);
% %(data4_0/max(data4_0))*-1
% title('4 Antennas with 20 degree');
% xlabel('Tx beam angle');ylabel('Normalized power');

%usrp 0度
% xxx = -90:10:90;
% data4_0_l = [64.221, 63.541, 63.041, 62.267, 61.798, 63.031, 63.977, 59.271, 58.147, 57.221];
% data4_0_r = [58.371, 59.246, 63.541, 63.107, 61.971, 62.414, 63.037, 63.841, 64.102];
% data4_0 =  [data4_0_l data4_0_r];
% subplot(1,2,2);
% plot(xxx,(data4_0/max(data4_0))*-1 +1.88);
% (data4_0/max(data4_0))*-1
% title('4 Antennas with 0 degree');
% xlabel('Tx beam angle');ylabel('Normalized power');

% %usrp 90度
% xxx = -90:10:90;
% data4_0_l = [65.725, 66.214, 67.968, 68.732, 69.867, 70.252, 71.537, 68.798, 68.322, 70.911];
% data4_0_r = [68.272, 68.564, 70.726, 68.726, 67.254, 66.308, 66.029, 65.926, 65.822];
% data4_0 =  [data4_0_l data4_0_r];
% subplot(1,2,2);
% plot(xxx,(data4_0/max(data4_0))*-1 +1.91);
% (data4_0/max(data4_0))*-1
% title('4 Antennas with 90 degree');
% xlabel('Tx beam angle');ylabel('Normalized power');

% %usrp 50度
% xxx = -90:10:90;
% data4_0_l = [62.046, 62.755, 64.542, 65.540, 66.876, 65.915, 63.962, 66.051, 65.841, 64.021];
% data4_0_r = [66.021, 66.241, 64.972, 63.211, 60.544, 61.721, 62.011, 62.431, 62.515];
% data4_0 =  [data4_0_l data4_0_r];
% subplot(1,2,2);
% plot(xxx,(data4_0/max(data4_0))*-1 +1.9);
% (data4_0/max(data4_0))*-1
% title('4 Antennas with 50 degree');
% xlabel('Tx beam angle');ylabel('Normalized power');
%%

