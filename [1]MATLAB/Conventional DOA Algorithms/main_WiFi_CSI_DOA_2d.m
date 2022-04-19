%% Programmed by ZZF: 2021.11.18
% This example is to practise the coding based on MUSIC algorithm

% clear;close all;clc
%% Tx antenna setting
nvar = 0.5; % noise power 0.01
method = 'frequency domain';    
cf = 2.45e9; % 10e6 2.45e9
lambda = physconst('LightSpeed') / cf;
M0 = 10; % Total number of antenna elements
nsnapshot = 800; % for wifi signals, only 80, 160, ... 800
K = 2; % # of source signal
angle = [-60,-30,20,50,50,70,-30,-45];
coherent_on = 1; % if make the sources coherent to each other
L = K; % # of subarrays for FSS
L_fb = K; % # of subarrays for FBSS   K/2
M = M0 - L + 1; % # of antenna elements in a subarray
m = M0-1; % # of antenna elements in ESPRIT subarrays
msg = ['The SNR is ',num2str(1/nvar), 'dB'];
disp(msg)

%% Help information
msg = ['The antenna freedom in FSS is ',num2str(M0-L+1-K)];
disp(msg)
msg = ['The antenna freedom in FBSS is ',num2str(M0-L_fb+1-K)];
disp(msg)
% if M0-L+1-K<1
%     msg = ['antenna number should be increased by at least ',...
%         num2str(1-(M0-L+1-K))];
%     disp(msg)
% end

msg = ['For spatial smoothing, ',num2str(2*K),...
    ' is the least required number of antenna element.'];
disp(msg)
msg = ['For forward backward spatial smoothing, ',num2str(3/2*K),...
    ' is the least required number of antenna element.'];
disp(msg)
% figure
% plot(s)
% title('Pulse');
% xlabel('Time (s)');
% ylabel('Amplitude (V)');

%% Configuret parameters of the WiFi signal and the channel 
% NOTE:If SNR is too low (lower than 20dB), the constelation 
% map is not very pretty!
flag_channel_on=1;flag_cfo_est=1;plot_on=0;
% rng(0) % Initialize the random number generator
%% Generate and decoder the CSI of the WiFi signal
snr = 30; %not used!
% %! JUST GET CSI FROM ONE OFDM SYMBOL
multiple_signal = zeros(56, K);
for i = 1: K
    [cfgVHT, tgacChannel] = Configurations(); %Configure the data field length to be 0.
    [tx]=generate_wifi_signal(cfgVHT,tgacChannel,snr,flag_channel_on);
    [IQsample,csi] = wifi_decoder(cfgVHT,tx,flag_cfo_est,plot_on,nsnapshot); 
    multiple_signal(:,i) = csi;
end
% if coherent_on
%     s = randn(nsnapshot,1)+1i*randn(nsnapshot,1);
%     s = [s,(1+2j)*s,(1-7j)*s,(4-2j)*s];
%     s = [s,(1+2j)*s];
% else
%     s = randn(nsnapshot,K)+1i*randn(nsnapshot,K);
% end
if isequal(method,'frequency domain')
    % NOTE: the length of the CSI is 56.
%     s = csi;
%     s = repmat(s, 1, K); % copy multiple input data for multiple sources
    s = multiple_signal;
else
    % NOTE: the lengtth of the IQ sample is 80, 
    % just the first OFDM symbol of the preamble field.
    s = IQsample; 
    s = repmat(s, 1, K); % copy multiple input data for multiple sources
end

%% Received signal at each array element
% Create the URA with default isotropic elements. 
% Set the frequency response range of the elements.
array = phased.URA('Size',[11 11],'ElementSpacing',[lambda/2 lambda/2]);
array.Element.FrequencyRange = [50.0e6 10.0e9];

doa1 = [-15;0];
doa2 = [17;25];

% Create the two signals and add random noise.
% t = (0:1/fs:1).';
% x1 = cos(2*pi*t*f1);
% x2 = cos(2*pi*t*f2);
x = collectPlaneWave(array,[s(:,1) s(:,2)],[doa1,doa2],fc);
noise = nvar/2*(randn(size(x))+1i*randn(size(x)));

% Create and execute the 2-D MUSIC estimator to find the directions of arrival.
estimator = phased.MUSICEstimator2D('SensorArray',array,...
    'OperatingFrequency',fc,...
    'NumSignalsSource','Property',...
    'ForwardBackwardAveraging', true,...
    'DOAOutputPort',true,'NumSignals',2,...
    'AzimuthScanAngles',-50:.5:50,...
    'ElevationScanAngles',-30:.5:30);
[~,doas] = estimator(x + noise)
% The estimated DOAs exactly match the true DOAs.

% Plot the 2-D spatial spectrum
plotSpectrum(estimator);


% % Add thermal noise (complex Gaussian distribution). \
% % The power of the noise is 0.5 watt. Here the power of the signal is 1 watt;
% % therefore, the SNR is 1/0.5 = 2, which is 3dB.
% rs = RandStream.create('mt19937ar','Seed',2008);
% noise = sqrt(nvar/2)*(randn(size(rx_ori)) + 1i*randn(size(rx_ori)));
% rx_noise = rx_ori + noise; % add noise to the signal
% 
% %% Find the spacial covariance matrix,Rxx, of the received signal
% % use the sample average hat{Rxx} to estimate the Rxx
% h_Rxx = rx_noise*rx_noise'/nsnapshot;
% 
% %% Various DOA estimation methods
% angle_scan = -90:0.01:90; % Scanning angle range
% [P_Bartlett]=Bartlett_DOA(h_Rxx,M0,lambda,angle_scan);
% % [P_LP]=LP_DOA(h_Rxx,M0,lambda,angle_scan);
% [P_ML]=ML_DOA(h_Rxx,M0,lambda,angle_scan);
% [P_capon]=Capon_DOA(h_Rxx,M0,lambda,angle_scan);
% [P_MUSIC]=MUSIC_DOA(h_Rxx,M0,lambda,angle_scan,K);
% [P_MUSIC_SS]=MUSIC_SS_DOA(h_Rxx,M0,L,lambda,angle_scan,K);
% [P_MUSIC_FBSS]=MUSIC_FBSS_DOA(h_Rxx,M0,L_fb,lambda,angle_scan,K);
% [est_angle]=ESPRIT_DOA(h_Rxx,m,K);
% % [est_angle]=ESPRIT_SS_DOA(h_Rxx,M0,m,K,L);
% 
% 
% %% Display the result
% figure
% plot(angle_scan,10*log10(P_MUSIC/max(P_MUSIC)),'Linewidth',1);
% % plot(angle_scan,10*log10(P_Bartlett/max(P_Bartlett)));
% grid on
% hold on
% % plot(angle_scan,10*log10(P_capon/max(P_capon)),'Linewidth',1);
% % plot(angle_scan,10*log10(P_LP/max(P_LP)),'Linewidth',1);
% % stem(est_angle,zeros(length(est_angle),1),'or','Linewidth',2);
% plot(angle_scan,10*log10(P_MUSIC_SS/max(P_MUSIC_SS)),'Linewidth',1);
% plot(angle_scan,10*log10(P_MUSIC_FBSS/max(P_MUSIC_FBSS)),'Linewidth',1);
% 
% for i = 1:K
%     msg = num2str(angle(i));
%     xline(angle(i),'--r',msg,'LabelOrientation','horizontal',...
%         'LabelHorizontalAlignment','center',...
%         'LabelVerticalAlignment','bottom')
% end
% 
% % legend('Bartlett','MVDR (Capon)','Linear Prediction','MUSIC','ESPRIT','MUSIC with SS','MUSIC with FBSS')
% legend('MUSIC','MUSIC with SS','MUSIC with FBSS')
% 
% title('DOA Estimation');
% xlabel('Angle (degree)');
% ylabel('Spatila Power Spectrum |P| (dB)');
