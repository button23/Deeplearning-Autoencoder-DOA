%%
%  To generate and prepare the DOA training data
%  Creator:  zzf  Date: 2022.01.07
%  Modified: zzf  Date: 2022.01.21
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
SNR = -20:20;
operSys = 'WINDOWS';
randn('seed',23);
s = sqrt(signalPower / 2) * randn(1, K) + sqrt(signalPower / 2) * 1i * randn(1, K);

%% Generate training and test datasets
% Create training data
[vec_noisy_train,vec_origin_train,Rxx_noisy_train,angle_doa_train] = DL_DataGenerator(s,'train',operSys,num_sample*5,SNR,M0,nsnapshot,K,lambda);
% Create test data
[vec_noisy_test,vec_origin_test,Rxx_noisy_test,angle_doa_test] = DL_DataGenerator(s,'test',operSys,num_sample,SNR,M0,nsnapshot,K,lambda);
%% Import the denoised data from the file
% Restore the SCM from denoised vector by autoencoder

% import the denoised vector
if isequal(operSys,'WINDOWS')
    dataPath = '/Users/button/OneDrive - 한양대학교/[7]Code/[1]Matlab/DOA_deeplearning/28-Jan-2022/data/result/denoised_data.mat'; % WINDOWS PATH:
else
    dataPath = '/Users/button/Deeplearning-Autoencoder-DOA/28-Jan-2022/data/result/denoised_data.mat'; % MAC PATH:
end
input_data = load(dataPath);
denoised_data = input_data.denoised_data;

% To restore real-valued vector to complex-valued sample covariance matrix.
[denoised_matrix] = DL_vec_2_SCM(denoised_data, num_sample, SNR, M0);
% isequal(denoised_matrix,Rxx_noisy_test)

%% Performance Comparison
angle_denoised = zeros(num_sample, K);
angle_noisy = zeros(num_sample, K);

angle_scan = -90:0.01:90; % Scanning angle range
% the performance of MUSIC using the denoised‰ scan

for i = 400:num_sample
    close all
    [P_MUSIC_denoised] = MUSIC_DOA(denoised_matrix(:, :,i), M0, lambda, angle_scan, K);
    [P_MUSIC_ori] = MUSIC_DOA(Rxx_noisy_test(:, :, i), M0, lambda, angle_scan, K);
    
    figure(1)
    plot(angle_scan, 10 * log10(P_MUSIC_denoised / max(P_MUSIC_denoised)), 'Linewidth', 1);
    grid on
    
    hold on
    plot(angle_scan, 10 * log10(P_MUSIC_ori / max(P_MUSIC_ori)), 'Linewidth', 1);
    msg = num2str(angle_doa_test(i, :));
    xline(angle_doa_test(i, :), '--r', msg, 'LabelOrientation', 'horizontal', ...
        'LabelHorizontalAlignment', 'right', ...
        'LabelVerticalAlignment', 'middle',...
        'LabelOrientation', 'horizontal',...
        'LineWidth', 1)
    %     msg = num2str(angle_doa_test(i,2));
    %     xline(angle_doa_test(i,2),'--r',msg,'LabelOrientation','horizontal',...
    %         'LabelHorizontalAlignment','center',...
    %         'LabelVerticalAlignment','middle')
    legend('denoised', 'noisy')
end

%%
figure
grid on
[~, locs] = findpeaks(P_MUSIC, angle_scan, 'Annotate', 'extents', 'Threshold', 0.0000001);

figure
grid on
plot(angle_scan, 10 * log10(P_MUSIC / max(P_MUSIC)), 'Linewidth', 1);
