% Generate the training/test data for autoencoder denosing DOA application
% DATE: 2022.01.07 by ZZF
%      2022.01.21 by ZZF
% Input:
%   signalPower:    The signal power
%   dataType:       Training data or test data
%   operSys         The operating system used by the computer
%   num_sample:     The number of samples
%   SNR:            The signal noise power ratio
%   M0:             Total number of antenna elements in the array
%   nsnapshot:      The number of snapshot
%   K:              The number of sources
%   lambda:         The wavelength
%   random_on:      Whether generate random sources for every sample
%
% Output:
%   vec_all:        The vector data vectorized from noisy SCM
%   vec_all_ori:    The vector data vectorized from noiseless SCM
%   Rxx_noisy       The noisy sample covariance matrix
%   angle_doa:      The random angle (# sample x # source)
%

% dada =load('denoised_data.mat','denoised_data')

function [vec_all,vec_all_ori,Rxx_noisy,angle_doa] = DL_DataGenerator(signalPower,dataType,operSys,num_sample,SNR,M0,nsnapshot,K,lambda, random_on)
%% The signal source
randn('seed',23);
s = sqrt(signalPower / 2) * randn(1, K) + sqrt(signalPower / 2) * 1i * randn(1, K);

%% Random Angle Generation
if K == 1
    ang1 = randi([-80, 80],num_sample,1);
    angle_doa = [ang1];
else
    ang1 = randi([-80, 5],num_sample,1);
    ang2 = randi([5, 80],num_sample,1);
    angle_doa = [ang1,ang2];
end
%% Generate the DOA Data
% h_Rxx = zeros(M0,M0); %% Here the bug ! this one needs to be moved to loop
vec_all = zeros(M0^2, num_sample);
vec_all_ori = zeros(M0^2, num_sample);
Rxx_noisy = zeros(M0,M0,num_sample);
for snr = 1 : length(SNR)
    tStart = tic;
    for a = 1: num_sample
        h_Rxx = zeros(M0,M0);
        nvar = signalPower/(10^(SNR(snr)/10)); % noise power
        %% Define incoming sources (multiply sqrt(2) to normalize the power)
        %% Received signal at each array element
        mu = -2*pi/lambda*lambda/2*sind(angle_doa(a,:));
        n = 0 : M0-1;
        % Steering matrix (considering the origin is in the center of the ULA)
        A = exp(1i*(n.'-(M0-1)/2)*mu);
        
        if random_on == 1 % controlled by the user
            s = sqrt(signalPower/ 2) * randn(1, K) + sqrt(signalPower / 2) * 1i * randn(1, K); % for random source
        end
        
        rx_ori = A*s.'; % Received signal at the array
        
        % Add thermal noise (complex Gaussian distribution)
        % The power of the noise is 0.5 watt. Here the power of the signal is 1 watt;
        % therefore, the SNR is 1/0.5 = 2, which is 3dB.
        %         rs = RandStream.create('mt19937ar','Seed',2008);
        for sn = 1:nsnapshot
            noise = sqrt(nvar/2)*(randn(size(rx_ori)) + 1i*randn(size(rx_ori)));
            rx_noise = rx_ori + noise; % add noise to the signal
            %% Find the spacial covariance matrix,Rxx, of the received signal
            % use the sample average hat{Rxx} to estimate the Rxx
            h_Rxx = h_Rxx + rx_noise*rx_noise'/nsnapshot;
        end
        
        Rxx = rx_ori*rx_ori'; % This is the SCM without noise.
        Rxx_noisy(:,:,a) = h_Rxx; % This is the SCM with noise.
        
        %% To convert complex-valued sample covariance matrix to real-valued vector
        [vec_all(:,a)]=DL_SCM_2_vec(h_Rxx,M0); % noisy SCM to vector
        [vec_all_ori(:,a)]=DL_SCM_2_vec(Rxx,M0); % noiseless SCM to vector
    end
    % Save to the files.
    DL_DataSave(dataType,operSys,SNR(snr),vec_all,vec_all_ori,angle_doa);
    
    tEnd = toc(tStart);
    fprintf('[%s] %d dB SNR Data out, consumed %d s\n',dataType,SNR(snr),tEnd);
end
end