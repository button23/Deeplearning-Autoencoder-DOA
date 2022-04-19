%% Programmed by ZZF: 2021.11.18
% This example is to practise the coding based on MUSIC algorithm

% clear;close all;clc
%% Tx antenna setting
nvar = 0.1; % noise power 0.01
method = 'tme domain';
cf = 2.45e9; % 10e6 2.45e9
lambda = physconst('LightSpeed') / cf;
M0 = 4; % Total number of antenna elements
nsnapshot = 80; % for wifi signals, only 80, 160, ... 800
K = 1; % # of source signal
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

%% get the DOA from the ray objects
load('rays_10cm.mat')
for i =1:3871
    for j = 1:4
        mid  = rays(j,i);
        
        
    end
end



%% Farse the data before applying the DOA estimation
single_AP = squeeze(features(:,:,1,:));
sigle_loc = squeeze(single_AP(:,:,1));
% sigle_loc = sigle_loc.';


%% Find the spacial covariance matrix,Rxx, of the received signal
% use the sample average hat{Rxx} to estimate the Rxx
h_Rxx = sigle_loc*sigle_loc'/24; % 24 samples for one CSI

%% Various DOA estimation methods
angle_scan = -90:0.01:90; % Scanning angle range
[P_Bartlett]=Bartlett_DOA(h_Rxx,M0,lambda,angle_scan);
% [P_LP]=LP_DOA(h_Rxx,M0,lambda,angle_scan);
[P_ML]=ML_DOA(h_Rxx,M0,lambda,angle_scan);
[P_capon]=Capon_DOA(h_Rxx,M0,lambda,angle_scan);
[P_MUSIC]=MUSIC_DOA(h_Rxx,M0,lambda,angle_scan,K);
[P_MUSIC_SS]=MUSIC_SS_DOA(h_Rxx,M0,L,lambda,angle_scan,K);
[P_MUSIC_FBSS]=MUSIC_FBSS_DOA(h_Rxx,M0,L_fb,lambda,angle_scan,K);
[est_angle]=ESPRIT_DOA(h_Rxx,m,K);
% [est_angle]=ESPRIT_SS_DOA(h_Rxx,M0,m,K,L);


%% Display the result
figure
plot(angle_scan,10*log10(P_MUSIC/max(P_MUSIC)),'Linewidth',1);
% plot(angle_scan,10*log10(P_Bartlett/max(P_Bartlett)));
grid on
hold on
% plot(angle_scan,10*log10(P_capon/max(P_capon)),'Linewidth',1);
% plot(angle_scan,10*log10(P_LP/max(P_LP)),'Linewidth',1);
% stem(est_angle,zeros(length(est_angle),1),'or','Linewidth',2);
plot(angle_scan,10*log10(P_MUSIC_SS/max(P_MUSIC_SS)),'Linewidth',1);
plot(angle_scan,10*log10(P_MUSIC_FBSS/max(P_MUSIC_FBSS)),'Linewidth',1);

for i = 1:K
    msg = num2str(angle(i));
    xline(angle(i),'--r',msg,'LabelOrientation','horizontal',...
        'LabelHorizontalAlignment','center',...
        'LabelVerticalAlignment','bottom')
end

% legend('Bartlett','MVDR (Capon)','Linear Prediction','MUSIC','ESPRIT','MUSIC with SS','MUSIC with FBSS')
legend('MUSIC','MUSIC with SS','MUSIC with FBSS')

title('DOA Estimation');
xlabel('Angle (degree)');
ylabel('Spatila Power Spectrum |P| (dB)');
