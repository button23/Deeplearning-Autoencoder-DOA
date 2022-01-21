% Generate the test data for autoencoder denosing DOA application
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

function [vec_all,vec_all_ori,Rxx_noisy,angle_doa]=DLtestDataGenerator(s,num_sample,SNR,M0,nsnapshot,K,lambda)
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
vec_all = zeros(M0^2, num_sample);
vec_all_ori = zeros(M0^2, num_sample);
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
        
        Rxx = rx_ori*rx_ori';
        
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
        
        Rxx_noisy(:,:,a) = h_Rxx;
    end
end
%% Prepare for the DL model
% Save Path
%     dataPath = '/home/hymc/bf_test/code_matlab/DOA_deeplearning/0107/test/test_data';
%     labelPath = '/home/hymc/bf_test/code_matlab/DOA_deeplearning/0107/test/test_label';
dataPath = '/Users/button/OneDrive - 한양대학교/[7]Code/[1]Matlab/DOA_deeplearning/0107/test/test_data';
originPath = '/Users/button/OneDrive - 한양대학교/[7]Code/[1]Matlab/DOA_deeplearning/0107/test/test_origin_data';
labelPath = '/Users/button/OneDrive - 한양대학교/[7]Code/[1]Matlab/DOA_deeplearning/0107/test/test_label';
%% Generate test data: Shuffle
test_data = vec_all.'; % make the data conform to numsample by (number antenna element^2(noisy))
% Generate random number
% rng(23)
% ind_ran = randperm(num_sample);
tarData = 'test_data';
save(dataPath, tarData)

%% Generate test data: Shuffle
test_origin_data = vec_all_ori.'; % make the data conform to numsample by (number antenna element^2(noisy))
tarData = 'test_origin_data';
save(originPath, tarData)

%% Generate label (DOA angle)
test_label = angle_doa;
tarData = 'test_label';
save(labelPath, tarData)

end
