% Compute the MUSIC Spatial Smoothing DOA estimate
% DATE: 2021.12.09 by ZZF
%
% Input:
%   h_Rxx:        The sample covariance matrix 
%   M0:           The total number of antenna elements in the array
%   L:            The number of subarrays
%   lambda:       The wavelength of the signal
%   angle_scan:   The scanning angles during angle search
%   K             The number of sources
% Output:
%   P_MUSIC_SS:   Spatial spectrum over scanning angles

function [P_MUSIC_SS]=MUSIC_SS_DOA(h_Rxx,M0,L,lambda,angle_scan,K)
%% Compute Spatial smoothing
mu = -2*pi/lambda*lambda/2*sind(angle_scan);
M = M0 - L + 1; % # of antenna elements in a subarray
ll = 0 : M-1;
% all steering vector M0 (# antenna element) x N (# of scan angle)
sv_all_ss = exp(1i*(ll.'-(M-1)/2)*mu);

% % The First way to calculate spatially smoothed covariance
% array_sub = 0;
% for i = 1 : L
%     rx_ori_sub = rx_noise(i:i+M-1,:);
%     array_sub = array_sub + rx_ori_sub*rx_ori_sub';
% end
% array_sub = array_sub /nsnapshot / L;

% %  SVD can directly work on the receivd signal matrix instead of covariance matrix
% [U,S,V] = svd(array_sub); 
%
% % Since there is only one source, only one large eigenvalue in S, and the
% % rest of the eigenvalues correspond to the noise eigenvectors.
% P = sum(S);                                     % power of received data
% P_cum = cumsum(S);                              % cumsum of V
% J = find(P_cum/P>=0.95);                        % or the coefficient is 0.9
% J = J(1);                                       % number of principal component
% V_n = U(:,J+1:end);

% The second way to calculate spatially smoothed covariance
% Find it directly from the calculated Rxx
array_sub_2 = 0;
for i = 1 : L
    array_sub_2 = array_sub_2 + h_Rxx(i:i+M-1,i:i+M-1);
end

% spatially smoothed covariance
array_sub_2 = array_sub_2 / L;

% Perform eigenvalue decomposition on Rxx to find the eigenvectors
% associated with the noise eigenvalues
[U,V] = eig(array_sub_2);
V = diag(V);             % vectorize eigenvalue matrix
[V,idx] = sort(V,'descend');                    % sort the eigenvalues in descending order
U = U(:,idx);                                   % reset the eigenvector
P = sum(V);                                     % power of received data
P_cum = cumsum(V);                              % cumsum of V

% define the noise space
J = find(P_cum/P>=0.95);                        % or the coefficient is 0.9
J = J(1);                                       % number of principal component
V_n = U(:,K+1:end);

P_MUSIC_SS = zeros(length(angle_scan),1);
for i = 1:length(angle_scan)
    sv = sv_all_ss(:,i);
    P_MUSIC_SS(i) = 1/abs(sv'* (V_n *V_n')*sv);
end

end