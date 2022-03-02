% Compute the MUSIC DOA estimate
% DATE: 2021.12.09 by ZZF
%
% Input:
%   h_Rxx:       The sample covariance matrix 
%   M0: the      Total number of antenna elements in the array
%   lambda:      The wavelength of the signal
%   angle_scan:  The scanning angles during angle search
%   K:          The number of sources
% Output:
%   P_MUSIC:     Spatial spectrum over scanning angles

function [P_MUSIC]=MUSIC_DOA(h_Rxx,M0,lambda,angle_scan,K)
mu = -2*pi/lambda*lambda/2*sind(angle_scan);
n = 0 : M0-1;
% all steering vector: M0 (# antenna element) x N (# of scan angle)
sv_all = exp(1i*(n.'-(M0-1)/2)*mu);

% Perform eigenvalue decomposition on Rxx to find the eigenvectors
% associated with the noise eigenvalues
[U,V] = eig(h_Rxx);
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
% [U,S,V] = svd(rx_noise);

% Since there is only one source, only one large eigenvalue in S, and the
% rest of the eigenvalues correspond to the noise eigenvectors.

% V_n = U(:,K+1:end);

P_MUSIC = zeros(length(angle_scan),1);
for i = 1:length(angle_scan)
    sv = sv_all(:,i);
    P_MUSIC(i) = 1/abs(sv'* (V_n *V_n')*sv);
end

end