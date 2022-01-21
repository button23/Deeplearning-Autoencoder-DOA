%%
%  To generate and prepare the DOA training data
%  Creator: zzf  Date: 2022.01.07
%
% clear
% clc
%% Tx antenna setting
cf = 2.45e9; % 10e6 2.45e9
lambda = physconst('LightSpeed') / cf;
M0 = 20; % Total number of antenna elements
nsnapshot = 1000;
K = 1; % # of source signal
num_sample = 1000; % number of samples
signalPower = 1; % Signal power
SNR = 10;

s = sqrt(signalPower/2)*randn(1,K)+sqrt(signalPower/2)*1i*randn(1,K);
%% Create training data
[~,angle_doa_train] = DLtrainDataGenerator(s,num_sample*5,SNR,M0,nsnapshot,K,lambda);

%% Create test data
[vec_noisy,vec_origin,Rxx_noisy_test,angle_doa_test] = DLtestDataGenerator(s,num_sample,SNR,M0,nsnapshot,K,lambda);

%% Performance Comparison
%% Restore the SCM from denoised vector by autoencoder
[denoised_matrix]=vec2SCM(num_sample,SNR,M0);

%% Performance Comparison
angle_denoised = zeros(num_sample,K);
angle_noisy = zeros(num_sample,K);

angle_scan = -90:0.01:90; % Scanning angle range
% the performance of MUSIC using the denoised‰ scan
for i = 1:num_sample
    close all
    [P_MUSIC_denoised]=MUSIC_DOA(denoised_matrix(:,:,i),M0,lambda,angle_scan,K);
    [P_MUSIC_ori]=MUSIC_DOA(Rxx_noisy_test(:,:,i),M0,lambda,angle_scan,K);

    figure(1)
    grid on
    plot(angle_scan,10*log10(P_MUSIC_denoised/max(P_MUSIC_denoised)),'Linewidth',1);
    hold on
    plot(angle_scan,10*log10(P_MUSIC_ori/max(P_MUSIC_ori)),'Linewidth',1);
    msg = num2str(angle_doa_test(i,1));
    xline(angle_doa_test(i,1),'--r',msg,'LabelOrientation','horizontal',...
        'LabelHorizontalAlignment','center',...
        'LabelVerticalAlignment','top')
%     msg = num2str(angle_doa_test(i,2));
%     xline(angle_doa_test(i,2),'--r',msg,'LabelOrientation','horizontal',...
%         'LabelHorizontalAlignment','center',...
%         'LabelVerticalAlignment','middle')
    legend('denoised','noisy')
end

%%
figure
grid on
[~,locs]=findpeaks(P_MUSIC,angle_scan,'Annotate','extents','Threshold',0.0000001);

figure
grid on
plot(angle_scan,10*log10(P_MUSIC/max(P_MUSIC)),'Linewidth',1);

