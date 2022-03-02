%%
%
%  Compare the performance of DOA between DL denoised sample covariance
%  matrix and the noisy sample covariance matrix
%  Creator: zzf  Date: 2022.01.07
%
% clear
% clc

targ = 'test';
%% Tx antenna setting
cf = 2.45e9; % 10e6 2.45e9
lambda = physconst('LightSpeed') / cf;
nsnapshot = 1000;
K = 2; % # of source signal

num_sample = 1000;
ang1 = randi([-80, 5],num_sample,1);
ang2 = randi([5, 80],num_sample,1);
angle = [ang1,ang2];

M0 = 20; % Total number of antenna elements
signalPower = 1; % Signal power

SNR = -10:10;
SNR = -20;

%% generate the test data
[Rxx_noisy,angle_doa]=testDataGenerator(num_sample,SNR,M0,nsnapshot,K,lambda);

%% Restore the SCM from denoised vector by autoencoder
[denoised_matrix]=vec2SCM(num_sample,SNR,M0);

%% Performance Comparison
angle_denoised = zeros(num_sample,K);
angle_noisy = zeros(num_sample,K);

angle_scan = -90:0.01:90; % Scanning angle range
% the performance of MUSIC using the denoised scan
for i = 1:num_sample
    close all
    [P_MUSIC_denoised]=MUSIC_DOA(denoised_matrix(:,:,i),M0,lambda,angle_scan,K);
    
    figure(1)
    grid on
    plot(angle_scan,10*log10(P_MUSIC_denoised/max(P_MUSIC_denoised)),'Linewidth',1);
    title('denoised')
    hold on
    msg = num2str(angle_doa(i,1));
    xline(angle_doa(i,1),'--r',msg,'LabelOrientation','horizontal',...
        'LabelHorizontalAlignment','center',...
        'LabelVerticalAlignment','top')
    msg = num2str(angle_doa(i,2));
    xline(angle_doa(i,2),'--r',msg,'LabelOrientation','horizontal',...
        'LabelHorizontalAlignment','center',...
        'LabelVerticalAlignment','middle')
    
    figure(2)
    [P_MUSIC_ori]=MUSIC_DOA(Rxx_noisy(:,:,i),M0,lambda,angle_scan,K);
    grid on
    plot(angle_scan,10*log10(P_MUSIC_ori/max(P_MUSIC_ori)),'Linewidth',1);
    title('noisy')
    hold on
    msg = num2str(angle_doa(i,1));
    xline(angle_doa(i,1),'--r',msg,'LabelOrientation','horizontal',...
        'LabelHorizontalAlignment','center',...
        'LabelVerticalAlignment','top')
    msg = num2str(angle_doa(i,2));
    xline(angle_doa(i,2),'--r',msg,'LabelOrientation','horizontal',...
        'LabelHorizontalAlignment','center',...
        'LabelVerticalAlignment','middle')
end

%%
figure
grid on
[~,locs]=findpeaks(P_MUSIC,angle_scan,'Annotate','extents','Threshold',0.0000001);

figure
grid on
plot(angle_scan,10*log10(P_MUSIC/max(P_MUSIC)),'Linewidth',1);
