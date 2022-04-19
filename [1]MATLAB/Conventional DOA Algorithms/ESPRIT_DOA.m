% Compute the ESPRIT DOA estimate
% DATE: 2021.12.09 by ZZF
%
% Input:
%   h_Rxx:      The sample covariance matrix 
%   m:          The total number of antenna elements in each subarray
%   K:          The number of signal sources

% Output:
%   est_angle:  Estimated angles

function [est_angle]=ESPRIT_DOA(h_Rxx,m,K)
% Perform eigenvalue decomposition on Rxx to find the eigenvectors
% associated with the noise eigenvalues
est_angle = 0;
[U,V] = eig(h_Rxx);
V = diag(V);             % vectorize eigenvalue matrix
[V,idx] = sort(V,'descend');                    % sort the eigenvalues in descending order
U = U(:,idx);                                   % reset the eigenvector
P = sum(V);                                     % power of received data
P_cum = cumsum(V);                              % cumsum of V
% define the signal subspace
indd= find(P_cum/P>=0.95);                        % or the coefficient is 0.9
P_n = indd(1);                                       % number of principal component
% if P_n == length(h_Rxx)
%    return 
% end
U_ss = U(:,1:K);                            % signal subspace eigenvectors

% Define selection matrix J
J_1 = [eye(m) zeros(m,1)];
J_2 = [zeros(m,1) eye(m)];
XX = J_1*U_ss;
YY = J_2*U_ss;

% invariance equation %J_1*U_ss*xxxx = J_2*U_ss;
Psi = estimation_TLS(XX,YY,K);

% Find the eigenvalues of Psi
Phi = eig(Psi);
arg = -1i*log(Phi);
theta_r = asin(-arg/pi);
est_angle = real(rad2deg(theta_r));
end