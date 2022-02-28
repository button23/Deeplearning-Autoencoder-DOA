%%
%  To generate and prepare the DOA training data
%  Creator:  zzf  Date: 2022.01.07
%  Modified: zzf  Date: 2022.01.21
%
% clear
% clc
%% Tx antenna setting
cf = 2.4e9; % 10e6 2.45e9
lambda = physconst('LightSpeed') / cf;
M0 = 8; % Total number of antenna elements # RF:test 8 # matlab 20
nsnapshot = 100;
K = 1; % # of source signal
num_sample = 500; % number of samples
signalPower = 1; % Signal power
SNR = -20:20;
operSys = 'UBUNTU';
random_on = 1;
%% Generate training and test datasets
% Create training data
[vec_noisy_train,vec_origin_train,Rxx_noisy_train,angle_doa_train] = DL_DataGenerator(signalPower,'train',operSys,num_sample*5,SNR,M0,nsnapshot,K,lambda,random_on);
% Create test data
[vec_noisy_test,vec_origin_test,Rxx_noisy_test,angle_doa_test] = DL_DataGenerator(signalPower,'test',operSys,num_sample,SNR,M0,nsnapshot,K,lambda,random_on);
%% Import the denoised data from the file
% Restore the SCM from denoised vector by autoencoder
basePath = '8_antenna_500_samples_100_snapshots';

% import the denoised vector
if isequal(operSys,'WINDOWS')
    dataPath = fullfile('C:\Users\HYPC300\OneDrive - 한양대학교\GitHub\Deeplearning-Autoencoder-DOA\data\result',basePath,'denoised_data.mat'); % WINDOWS PATH:
elseif isequal(operSys,'UBUNTU')
    dataPath = fullfile('/home/hymc/[0]Github/Deeplearning-Autoencoder-DOA/[2]Tensorflow/[1]DAE_DOA_implementation/result',basePath,'denoised_data.mat'); % UBUNTU PATH: 
else
    dataPath = '/Users/button/Deeplearning-Autoencoder-DOA/28-Jan-2022/data/result/denoised_data.mat'; % MAC PATH:
end
input_data = load(dataPath);
denoised_data = input_data.denoised_data;

% To restore real-valued vector to complex-valued sample covariance matrix.
% [denoised_matrix] = DL_vec_2_SCM(denoised_data, num_sample, SNR, M0);
[denoised_matrix] = DL_vec_2_SCM(denoised_data, 1, SNR, M0); % delete

% isequal(denoised_matrix,Rxx_noisy_test)

%% Import the test data and label (DOA angle)
testPath = fullfile('/home/hymc/[0]Github/Deeplearning-Autoencoder-DOA/[2]Tensorflow/[1]DAE_DOA_implementation/result',basePath,'test_data.mat'); % WINDOWS PATH:
input_data = load(testPath);
test_data = input_data.test_data;

% [Rxx_noisy_test] = DL_vec_2_SCM(test_data, num_sample, SNR, M0);
[Rxx_noisy_test] = DL_vec_2_SCM(test_data, 1, SNR, M0); % to delete


%%
labelPath = fullfile('/home/hymc/[0]Github/Deeplearning-Autoencoder-DOA/[2]Tensorflow/[1]DAE_DOA_implementation/result',basePath,'test_label.mat'); % WINDOWS PATH:
input_data = load(labelPath);
angle_doa_test = input_data.test_label;

%% Performance Comparison
angle_denoised = zeros(num_sample, K);
angle_noisy = zeros(num_sample, K);

angle_scan = -90:0.01:90; % Scanning angle range
% the performance of MUSIC using the denoised‰ scan

for i = 1:1000
    close all
    [P_MUSIC_denoised] = MUSIC_DOA(denoised_matrix(:, :,i), M0, lambda, angle_scan, K);
    [P_MUSIC_ori] = MUSIC_DOA(Rxx_noisy_test(:, :, i), M0, lambda, angle_scan, K);

%     [P_MUSIC_denoised] = MUSIC_DOA(denoised_matrix, M0, lambda, angle_scan, K); % delete
%     [P_MUSIC_ori] = MUSIC_DOA(h_Rxx, M0, lambda, angle_scan, K); % delete

    figure(1)
    plot(angle_scan, 10 * log10(P_MUSIC_denoised / max(P_MUSIC_denoised)), 'Linewidth', 1);
    grid on
    
    hold on
    plot(angle_scan, 10 * log10(P_MUSIC_ori / max(P_MUSIC_ori)), 'Linewidth', 1);
       msg = num2str(0); % delete
    xline(0, '--r', msg, 'LabelOrientation', 'horizontal', ...
        'LabelHorizontalAlignment', 'right', ...
        'LabelVerticalAlignment', 'middle',...
        'LabelOrientation', 'horizontal',...
        'LineWidth', 1)  % delete
    
%     msg = num2str(angle_doa_test(i,1));
%     xline(angle_doa_test(i, 1), '--r', msg, 'LabelOrientation', 'horizontal', ...
%         'LabelHorizontalAlignment', 'right', ...
%         'LabelVerticalAlignment', 'middle',...
%         'LabelOrientation', 'horizontal',...
%         'LineWidth', 1)
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
