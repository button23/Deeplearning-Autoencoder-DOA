% Compute the ESPRIT Spatial Spacing DOA estimate
% DATE: 2021.12.09 by ZZF
%
% Input:
%   h_Rxx:      The sample covariance matrix 
%   M0:         The total number of antenna elements in array
%   m:          The total number of antenna elements in each subarray
%   K:          The number of signal sources
%   L:          The number of subarrays for spatial smoothing

% Output:
%   est_angle:  Estimated angles

function [est_angle_ss]=ESPRIT_SS_DOA(h_Rxx,M0,m,K,L)
% The second way to calculate spatially smoothed covariance
% Find it directly from the calculated Rxx
M = M0 - L + 1; % # of antenna elements in a subarray

array_fb = 0;
for i = 1 : L
    array_fb = array_fb + h_Rxx(i:i+M-1,i:i+M-1);
end

% spatially smoothed covariance
array_fb = array_fb / L;

% Perform eigenvalue decomposition on Rxx to find the eigenvectors
% associated with the noise eigenvalues
[U,V] = eig(array_fb);
V = diag(V);             % vectorize eigenvalue matrix
[V,idx] = sort(V,'descend');                    % sort the eigenvalues in descending order
U = U(:,idx);                                   % reset the eigenvector
P = sum(V);                                     % power of received data
P_cum = cumsum(V);                              % cumsum of V
% define the signal subspace
indd= find(P_cum/P>=0.95);                        % or the coefficient is 0.9
P_n = indd(1);                                       % number of principal component
U_ss = U(:,1:P_n);                            % signal subspace eigenvectors

% Define selection matrix J
J_1 = [eye(m) zeros(m,1)];
J_2 = [zeros(m,1) eye(m)];

% invariance equation %J_1*U_ss*xxxx = J_2*U_ss;

% TLS
XX = J_1*U_ss;
YY = J_2*U_ss;
z = [XX YY]; % put matrix x and vector y together
[~ , S_z, V_z]=svd(z);
n = length(XX(1,:)); % number o columns
Psi = -V_z(1:n,n+1:end)*inv(V_z(n+1:end,n+1:end)); % the most important line!!!
% xtyt = -z*V(:,n+1)*V(:,n+1)'; % Find the error of x and y (x_tilda and y_tilda)
%
% % Real Line fitting with the estimated parameter a_tls
% y_tls2 = (XX+xtyt(:,1:n))*a_tls; % (x+x_tilda) multiply by the parameter a_tls

% Find the eigenvalues of Psi
Phi = eig(Psi);
arg = -1i*log(Phi);
theta_r = asin(-arg/pi);
est_angle_ss = real(rad2deg(theta_r));
end