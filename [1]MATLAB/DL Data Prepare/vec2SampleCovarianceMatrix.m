% Restore the sample covariance matrix based on the denoised vector
% by autoencoder.
% DATE: 2022.01.07 by ZZF
%
% Input:
%   num_sample:         The number of samples
%   SNR:                The signal noise power ratio
%   M0:                 Total number of antenna elements in the array
%   nsnapshot:          The number of snapshot
%   K:                  The number of sources
% Output:
%   denoised_matrix     The denoised sample covariance matrix


function [denoised_matrix]=vec2SCM(num_sample,SNR,M0)
    % import the denoised vector 
    dataPath = '/Users/button/OneDrive - 한양대학교/[7]Code/[1]Matlab/DOA_deeplearning/0107/result/denoised_data.mat';
    input_data = load(dataPath);
    denoised_data = input_data.denoised_data;
   
    %% DL denoised sample covariance matrix based music
    denoised_matrix = zeros(M0,M0,num_sample);
    for snr = 1 : length(SNR)
        for a = 1: num_sample
            %restore the real part
            oneshot = denoised_data(a,:);
            ones_low_real = tril(ones(M0,M0));
            index_low = logical(ones_low_real); % change to logical value for indexing
            ones_low_real(index_low)=oneshot(1:M0*(M0+1)/2);

            % restoure the imaginary part
            ones_low_imag = tril(ones(M0,M0),-1);
            index_low = logical(ones_low_imag); % change to logical value for indexing
            ones_low_imag(index_low)=oneshot(M0*(M0+1)/2+1:end);

            % combine real part and imaginary part
            low_part = ones_low_real + 1i*ones_low_imag;

            % to generate the high part
            high_part = tril(low_part,-1)';

            % combine the low and high triangular part
            denoised_matrix(:,:,a) = low_part + high_part;
        end
    end
            denoised_matrix = conj(denoised_matrix);
end