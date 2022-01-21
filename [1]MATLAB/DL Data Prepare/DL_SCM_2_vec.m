% Convert complex-valued sample covariance matrix to real-valued vector
% DATE: 2022.01.21 by ZZF
%
% Input:
%   Rxx:           The noiseless Sample Covariance Matrix
%   h_Rxx:         The noisy Sample Covariance Matrix
%   num_sample:    The number of samples
%   M0:            Total number of antenna elements in the array
%
% Output:
%   Rxx_noisy      The noisy sample covariance matrix
%   angle_doa:     The random angle (# sample x # source)
%   

function [vec_all,vec_all_ori]=DL_SCM_2_vec(Rxx,h_Rxx,num_sample,M0)
        vec_all = zeros(M0^2, num_sample);
        vec_all_ori = zeros(M0^2, num_sample);
        %% To convert complex-valued sample covariance matrix to real-valued vector
        h_Rxx_triu = triu(h_Rxx); % get the up-triangular matrix of noised SCM
        Rxx_triu = triu(Rxx); % get the up-triangular matrix of pure SCM
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