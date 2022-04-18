%%
%  To generate and prepare the DOA training data
%  Creator: zzf  Date: 2022.01.07
%
% clear
% clc
targ = 'test';
%% Tx antenna setting
cf = 2.45e9; % 10e6 2.45e9
lambda = physconst('LightSpeed') / cf;
M0 = 20; % Total number of antenna elements
nsnapshot = 1000;
K = 2; % # of source signal

num_sample = 1000;
ang1 = randi([-80, 5],num_sample,1);
ang2 = randi([5, 80],num_sample,1);
angle = [ang1,ang2];

coherent_on = 0; % if make the sources coherent to each other
L = K; % # of subarrays for FSS
L_fb = K; % # of subarrays for FBSS   K/2
M = M0 - L + 1; % # of antenna elements in a subarray
m = M0-1; % # of antenna elements in ESPRIT subarrays
signalPower = 1; % Signal power

SNR = -10:10;
SNR = -20;

%% Generate the DOA Data
vec_all = zeros(M0^2, num_sample);
vec_all_ori = zeros(M0^2, num_sample);

for snr = 1 : length(SNR)
    for a = 1: num_sample
        nvar = signalPower/(10^(SNR(snr)/10)); % noise power
        %% Define incoming sources (multiply sqrt(2) to normalize the power)
        s = sqrt(signalPower/2)*randn(nsnapshot,K)+sqrt(signalPower/2)*randn(nsnapshot,K);
        %% Received signal at each array element
        mu = -2*pi/lambda*lambda/2*sind(angle(a,:));
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
        Rxx = rx_ori*rx_ori'/nsnapshot;
        %% To convert complex-valued sample covariance matrix to real-valued vector
        h_Rxx_triu = triu(h_Rxx);
        Rxx_triu = triu(Rxx);
        h_Rxx_tril = h_Rxx_triu.'; % convert to lower triangular matirx (in order to extract value in a desired manner)
        Rxx_tril = Rxx_triu.';
        
        ones_low = tril(ones(M0,M0));
        ones_low = logical(ones_low); % change to logical value for indexing
        vec_rxx = h_Rxx_tril(ones_low); % extract values from the corresponding index for vectorization
        vec_rxx_ori = Rxx_tril(ones_low);
        vec_real = real(vec_rxx); % extract the real part of the vector (N(N+1)/2 entries)
        vec_real_ori = real(vec_rxx_ori);
        
        ones_low = tril(ones(M0,M0),-1); % do not include the diagonal
        ones_low = logical(ones_low); % change to logical value for indexing
        vec_rxx = h_Rxx_tril(ones_low); % extract values from the corresponding index for vectorization
        vec_rxx_ori = Rxx_tril(ones_low);
        
        vec_imag = imag(vec_rxx); % extract the imag part of the vector (N(N-1)/2 entries)
        vec_imag_ori = imag(vec_rxx_ori);
        vec_one = [vec_real; vec_imag]; % combine the real and imaginary part (N * N entries)
        vec_one_ori = [vec_real_ori; vec_imag_ori];
        
        vec_all(:,a) = vec_one;
        vec_all_ori(:,a) = vec_one_ori;
    end
end
%% Prepare for the DL model
% Load Path and Save Path
tarData = strcat(targ, '_data');
tarLabel = strcat(targ, '_label');

dataPath = strcat(pwd,'/0107/');
dataPath = strcat(dataPath, targ, '/', tarData);
labelPath = strcat(pwd,'/0107/');
labelPath = strcat(labelPath, targ, '/', tarLabel);

%% Generate training data: Shuffle
rx_noise_ = vec_all.'; % make the data conform to numsample by (number antenna element^2(noisy))
% Generate random number
% rng(23)
% ind_ran = randperm(num_sample);

if isequal(targ, 'train')
    train_data = rx_noise_;
else
    test_data = rx_noise_;
end
save(dataPath, tarData)

%% Generate label (original data)
rx_origin = vec_all_ori.'; % make the data conform to numsample by (number antenna element^2(original))

if isequal(targ, 'train')
    train_label = rx_origin;
else
    test_label = rx_origin;
end
save(labelPath, tarLabel)
