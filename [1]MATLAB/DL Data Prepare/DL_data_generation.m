%%
%
%  To generate and prepare the DOA training data
%  Creator: zzf  Date: 2021.12.24
%
% clear
% clc
targ = 'train';
%% Tx antenna setting
cf = 2.45e9; % 10e6 2.45e9
lambda = physconst('LightSpeed') / cf;
M0 = 20; % Total number of antenna elements
nsnapshot = 1000;
K = 2; % # of source signal
angle1 = -80:-5;
angle2 = 5:80;
coherent_on = 0; % if make the sources coherent to each other
L = K; % # of subarrays for FSS
L_fb = K; % # of subarrays for FBSS   K/2
M = M0 - L + 1; % # of antenna elements in a subarray
m = M0-1; % # of antenna elements in ESPRIT subarrays
signalPower = 1; % Signal power

SNR = -10:10;
SNR = -20;

num_class = length(angle);
%% Generate the DOA Data
rx_all = zeros(M0, nsnapshot*num_class);
for snr = 1 : length(SNR)
    for a = 1: num_class
        nvar = signalPower/(10^SNR(snr)); % noise power
        %% Define incoming sources (multiply sqrt(2) to normalize the power)
        s = sqrt(signalPower/2)*randn(nsnapshot,1)+sqrt(signalPower/2)*randn(nsnapshot,1);
        %% Received signal at each array element
        mu = -2*pi/lambda*lambda/2*sind(angle(a(1:K)));
        n = 0 : M0-1;
        % Steering matrix (considering the origin is in the center of the ULA)
        A = exp(1i*(n.'-(M0-1)/2)*mu);
        
        rx_ori = A*s.'; % Received signal at the array
        
        % Add thermal noise (complex Gaussian distribution)
        % The power of the noise is 0.5 watt. Here the power of the signal is 1 watt;
        % therefore, the SNR is 1/0.5 = 2, which is 3dB.
%         rs = RandStream.create('mt19937ar','Seed',2008);
        
        noise = sqrt(nvar/2)*(randn(size(rx_ori)) + 1i*randn(size(rx_ori)));
        rx_noise = rx_ori + noise; % add noise to the signal
        %% Find the spacial covariance matrix,Rxx, of the received signal
        % use the sample average hat{Rxx} to estimate the Rxx
        h_Rxx = rx_noise*rx_noise'/nsnapshot;
        rx_all(:,(a-1)*nsnapshot+1:a*nsnapshot) = rx_noise;
    end
end
%% Prepare for the DL model
% Load Path and Save Path

tarData = strcat(targ, '_data');
tarLabel = strcat(targ, '_label');

dataPath = strcat(pwd,'/1224/');
dataPath = strcat(dataPath, targ, '/', tarData);
labelPath = strcat(pwd,'/1224/');
labelPath = strcat(labelPath, targ, '/', tarLabel);

%% Generate training data: Shuffle
rx_noise_ = rx_all.'; % make the data conform to numsample by numant_ele
% Generate random number
rx_noise_real = real(rx_noise_);
rx_noise_imag = imag(rx_noise_);
rx_noise_split = [rx_noise_real, rx_noise_imag];

rng(23)
ind_ran = randperm(nsnapshot * num_class);
if isequal(targ, 'train')
    train_data = rx_noise_split(ind_ran',:,:);
else
    test_data = rx_noise_split(ind_ran',:,:);
end
save(dataPath, tarData)

%% Generate Lable
label = zeros(nsnapshot, 1);
ss = 0;

for k = 1 : num_class
    label((k-1)*nsnapshot+1:k*nsnapshot,:) = ones(nsnapshot,1) * ss;
    ss = ss + 1;
end
if isequal(targ, 'train')
    train_label = label(ind_ran);
else
    test_label = label(ind_ran);
end
save(labelPath, tarLabel)

