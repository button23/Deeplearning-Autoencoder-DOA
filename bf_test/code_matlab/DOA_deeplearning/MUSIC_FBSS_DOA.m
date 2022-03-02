% Compute the MUSIC Forward Backward Spatial Smoothing DOA estimate
% DATE: 2021.12.09 by ZZF
%
% Input:
%   h_Rxx:          The sample covariance matrix 
%   M0:             The total number of antenna elements in the array
%   L_fb:           The number of subarrays
%   lambda:         The wavelength of the signal
%   angle_scan:     The scanning angles during angle search
%   K               The number of sources
% Output:
%   P_MUSIC_FBSS:   Spatial spectrum over scanning angles

function [P_MUSIC_FBSS]=MUSIC_FBSS_DOA(h_Rxx,M0,L_fb,lambda,angle_scan,K)
%% Compute Spatial smoothing
mu = -2*pi/lambda*lambda/2*sind(angle_scan);
M = M0 - L_fb + 1; % # of antenna elements in a subarray
ll = 0 : M-1;
% All steering vector: M0 (# antenna element) x N (# of scan angle)
sv_all_ss = exp(1i*(ll.'-(M-1)/2)*mu);

% % Fisrt method to just compute the backward spatially smoothed covariance
% % matrix
% array_sub_back = 0;
% for i = 1 : L_fb
%     rx_noise_ud = flipud(rx_noise);
%     rx_ori_sub = conj(rx_noise_ud(i:i+M-1,:));
%     array_sub_back = array_sub_back + rx_ori_sub*rx_ori_sub';
% end
% %
% array_sub_back = array_sub_back / L_fb/nsnapshot;
%
% % % the forward/backward smoothed covariance matrix
% R_fb = (array_sub + array_sub_back)/2;


% Second method to compute the forward backward spatially smoothed
% covariance matrix at once
H = fliplr(eye(M0));
% Select subarray elements directly from Rxx after updown and leftright inverse
h_Rxxb = H*conj(h_Rxx)*H;
array_sub_fb = 0;
h_Rxxfb = (h_Rxx + h_Rxxb)/2;

for i = 1 : L_fb
    array_sub_fb = array_sub_fb + h_Rxxfb(i:i+M-1,i:i+M-1);
end

% spatially smoothed covariance
array_sub_fb = array_sub_fb / L_fb;


% Compute the MUSIC DOA estimate
% Perform eigenvalue decomposition on Rxx to find the eigenvectors
% associated with the noise eigenvalues
[U,V] = eig(array_sub_fb);
V = diag(V);             % vectorize eigenvalue matrix
[V,idx] = sort(V,'descend');                    % sort the eigenvalues in descending order
U = U(:,idx);                                   % reset the eigenvector
P = sum(V);                                     % power of received data
P_cum = cumsum(V);                              % cumsum of V

% define the noise space
J = find(P_cum/P>=0.95);                        % or the coefficient is 0.9
J = J(1);                                       % number of principal component
V_n = U(:,K+1:end);

% % SVD can directly work on the receivd signal matrix instead of covariance matrix
% [U,S,V] = svd(array_sub_fb);

% Since there is only one source, only one large eigenvalue in S, and the
% rest of the eigenvalues correspond to the noise eigenvectors.
% V_n = U(:,K+1:end);

P_MUSIC_FBSS = zeros(length(angle_scan),1);
for i = 1:length(angle_scan)
    sv = sv_all_ss(:,i);
    P_MUSIC_FBSS(i) = 1/abs(sv'* (V_n *V_n')*sv);
end

end