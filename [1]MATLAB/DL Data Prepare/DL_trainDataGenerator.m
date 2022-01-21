% Generate the training data for autoencoder denosing DOA application
% DATE: 2022.01.07 by ZZF
%
% Input:
%   num_sample:      The number of samples
%   SNR:          The signal noise power ratio
%   M0:           Total number of antenna elements in the array
%   nsnapshot:    The number of snapshot
%   K:            The number of sources
%   lambda:       The wavelength
% Output:
%   Rxx_noisy     The noisy sample covariance matrix
%   angle_doa:        The random angle (# sample x # source)
%

% dada =load('denoised_data.mat','denoised_data')

function [Rxx_noisy,angle_doa]=DL_trainDataGenerator(s,num_sample,SNR,M0,nsnapshot,K,lambda)
signalPower = 1;
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
h_Rxx = zeros(M0,M0);

Rxx_noisy = zeros(M0,M0,num_sample);
for snr = 1 : length(SNR)
    for a = 1: num_sample
        nvar = signalPower/(10^(SNR(snr)/10)); % noise power
        %% Define incoming sources (multiply sqrt(2) to normalize the power)
        %% Received signal at each array element
        mu = -2*pi/lambda*lambda/2*sind(angle_doa(a,:));
        n = 0 : M0-1;
        % Steering matrix (considering the origin is in the center of the ULA)
        A = exp(1i*(n.'-(M0-1)/2)*mu);
        
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
        Rxx_noisy(:,:,a) = h_Rxx;
        
        %% To convert complex-valued sample covariance matrix to real-valued vector
        [vec_all,vec_all_ori]=DL_SCM_2_vec(Rxx,h_Rxx,num_sample,M0);
        
    end
end
%% Prepare for the DL model
dataPath = 'C:\Users\HYPC300\OneDrive - 한양대학교\GitHub\Deeplearning-Autoencoder-DOA\data\2022_01_21\train\train_data';
originPath = 'C:\Users\HYPC300\OneDrive - 한양대학교\GitHub\Deeplearning-Autoencoder-DOA\data\2022_01_21\train\origin_data';
labelPath = 'C:\Users\HYPC300\OneDrive - 한양대학교\GitHub\Deeplearning-Autoencoder-DOA\data\2022_01_21\train\train_label';
%% Generate test data: Shuffle
train_data = vec_all.'; % make the data conform to numsample by (number antenna element^2(noisy))
tarData = 'train_data';
save(dataPath, tarData)

origin_data = vec_all_ori.'; % make the data conform to numsample by (number antenna element^2(noisy))
tarData = 'origin_data';
save(originPath, tarData)

%% Generate label (DOA angle)
train_label = angle_doa;
tarData = 'train_label';
save(labelPath, tarData)