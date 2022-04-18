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
%%
% 
txarraystv = phased.SteeringVector('SensorArray',txarray,'PropagationSpeed',c); 
txang = [60;0];  % tx angle (azimuth; elevation)          -180<=azimuth<=180  -90<=elevation<=90
wt = txarraystv(fc,txang)';
% wt = wt.'


txbeam_ang = -90:90;
%txbeam_ang = -180:180;
txbeam_ang_rad = (pi*txbeam_ang)/180;
txsv = steervec(txmipos(2,:),txbeam_ang);
txbeam = abs(wt*txsv);  
txbeam(91)
%txbeam = txbeam/max(txbeam);
%[txbeampos_x,txbeampos_y] = pol2cart(deg2rad(txbeam_ang),txbeam);

% hFig = figure(1);
% set(hFig, 'Position', [0 0 800 400]);
figure()
% subplot(1,2,1);
plot(txbeam_ang,txbeam);
xlabel('Tx Beam Angle');ylabel('Normalized Amplitude');
set(gca,'xtick',-180:10:180)
%xlim([txbeam_ang(1) txbeam_ang(end)]); ylim([0 1.0]);
title('Simulation result')

% wt = wt.'
%% weight multiplication
data_raw = var.waveform(1:15360);
data_raw_1 = data_raw.*wt(1);
data_raw_2 = data_raw.*wt(2);
data = [data_raw_1 data_raw_2];
data_scale = data * 2^15;
data_scale_short = int16(data_scale);
data_scale_short_real = real(data_scale_short);
data_scale_short_imag = imag(data_scale_short);

data_real1 = data_scale_short_real(:,1).';
data_real2 = data_scale_short_real(:,2).';
data_imag1 = data_scale_short_imag(:,1).';
data_imag2 = data_scale_short_imag(:,2).';

data1 = [data_imag1;data_real1];
data2 = [data_imag2;data_real2];

data1_res = reshape(data1,30720,1);
data2_res = reshape(data2,30720,1);

sent_data = [data1_res,data2_res];

%% usrp setting
%%
%   Programmer: ZZF
%   Date: 2021.9.10
%
%% Tx USRP Initialization
% release(radio);
radio = comm.SDRuTransmitter(...
    'Platform',             'X300', ...
    'IPAddress',            '192.168.1.23', ...
    'MasterClockRate',      184.32e6, ...
    'CenterFrequency',      2.4e9, ...
    'Gain',                 20, ...
    'InterpolationFactor',  1, ...
    'ClockSource',          'External', ...
    'PPSSource',            'External', ...
    'ChannelMapping',       [1 2]);




%%
% Inside a while loop, transmit the sine wave using the tx System object.
% Display a message when transmission is complete. Release the radio System object.
while true
    radio(sent_data);
end
%
% disp("Transmission Stopped");

