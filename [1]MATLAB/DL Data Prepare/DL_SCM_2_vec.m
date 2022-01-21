% Convert complex-valued sample covariance matrix to real-valued vector
% DATE: 2022.01.21 by ZZF
%
% Input:
%   Rxx:           The Sample Covariance Matrix
%   M0:            Total number of antenna elements in the array
%
% Output:
%   vec_one        The vector converted from SCM

function [vec_one]=DL_SCM_2_vec(Rxx,M0)

%% To convert complex-valued sample covariance matrix to real-valued vector
Rxx_triu = triu(Rxx); % get the up-triangular matrix of pure SCM
Rxx_tril = Rxx_triu.';% convert to lower triangular matirx (in order to extract value in a desired manner)

% Get the real part of half of the all the element in SCM (including the diagonal)
ones_low = logical(tril(ones(M0,M0))); % change to logical value for indexing
vec_rxx = Rxx_tril(ones_low); % extract values from the corresponding index (not including diagnal) for vectorization
vec_real = real(vec_rxx); % extract the real part of the vector (M0(M0+1)/2 entries)

% Get the imaginary part of half of the all the element in SCM (not including the diagonal)
ones_sublow = logical(tril(ones(M0,M0),-1)); % change to logical value for indexing (% do not include the diagonal)
vec_rxx = Rxx_tril(ones_sublow); % extract values from the corresponding index for vectorization
vec_imag = imag(vec_rxx); % extract the imag part of the vector (M0(M0-1)/2 entries)

vec_one = [vec_real; vec_imag]; % combine the real and imaginary part (N * N entries)

end