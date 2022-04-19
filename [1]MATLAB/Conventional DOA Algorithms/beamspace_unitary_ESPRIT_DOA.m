% Compute the Beamspace Unitary ESPRIT DOA estimate
% DATE: 2021.12.15 by ZZF
%
% Input:
%   h_Rxx:      The sample covariance matrix
%   m:          The total number of antenna elements in each subarray
%   M0:          The total number of antenna elements in the array
%   K:          The number of signal sources

% Output:
%   est_angle:  Estimated angles

function [est_angle]=beamspace_unitary_ESPRIT_DOA(h_Rxx,m,M0,K)
%% DFT Beamspace
% selection matrix Gamma1
gamma1 = zeros()

% selection matrix Gamma1


%% Second method to compute the forward backward spatially smoothed
% covariance matrix at once
H = fliplr(eye(M0));
% Select subarray elements directly from Rxx after updown and leftright inverse
h_Rxxb = H*conj(h_Rxx)*H;
h_Rxxfb = h_Rxx + h_Rxxb;

% Define left Pi real symmetric matrix Q_m
if mod(m,2)==0 % if m is even
    Q_m = 1/sqrt(2) * [eye(m/2) 1i*eye(m/2); fliplr(eye(m/2)) -1i*fliplr(eye(m/2))];
else % else if m is odd
    mm = m-1;
    Q_m = 1/sqrt(2) * [eye(mm/2) zeros(mm/2,1) 1i*eye(mm/2); zeros(1,mm/2) sqrt(2) zeros(1,mm/2);fliplr(eye(mm/2)) zeros(mm/2,1) -1i*fliplr(eye(mm/2))];
end
% Define left Pi real symmetric matrix Q_M
if mod(M0,2)==0 % if M0 is even
    Q_M = 1/sqrt(2) * [eye(M0/2) 1i*eye(M0/2); fliplr(eye(M0/2)) -1i*fliplr(eye(M0/2))];
else % else if M0 is odd
    MM0 = M0-1;
    Q_M = 1/sqrt(2) * [eye(MM0/2) zeros(MM0/2,1) 1i*eye(MM0/2); zeros(1,MM0/2) sqrt(2) zeros(1,MM0/2);fliplr(eye(MM0/2)) zeros(MM0/2,1) -1i*fliplr(eye(MM0/2))];
end

% Calculate the real-valued subspace estimation E_s
GammaX = Q_M'*h_Rxxfb*Q_M;

%% SVD method using the received data (Direct)
% % Define left Pi real symmetric matrix Q_M
% Q_2N = 1/sqrt(2) * [eye(N) 1i*eye(N); fliplr(eye(N)) -1i*fliplr(eye(N))];
% size_mat = size(h_Rxx);
% H1 = fliplr(eye(size_mat(1)));
% H2 = fliplr(eye(size_mat(2)));
%
% Z = [h_Rxx H1*h_Rxx*H2];
% GammaX = Q_M'*Z*Q_2N;
%
%
% [U,S,V] = svd(GammaX);
% V_n = U(:,K+1:end);
%% EVD method using Covariance matrix
% Perform eigenvalue decomposition on Rxx to find the eigenvectors
% associated with the noise eigenvalues
[U,V] = eig(GammaX);
V = diag(V);                                    % vectorize eigenvalue matrix
[V,idx] = sort(V,'descend');                    % sort the eigenvalues in descending order
U = U(:,idx);                                   % reset the eigenvector
P = sum(V);                                     % power of received data
P_cum = cumsum(V);                              % cumsum of V
% define the signal subspace
indd= find(P_cum/P>=0.95);                      % or the coefficient is 0.9
P_n = indd(1);                                  % number of principal component
E_s = U(:,1:K);       % P_n K                     % real-valued signal subspace eigenvectors

% Define selection matrix J
J_1 = [eye(m) zeros(m,1)];
J_2 = [zeros(m,1) eye(m)];

% K1 (real part) and K2 (imaginary part)
% K1 = 1/2*Q_m'*(J_1+J_2)*Q_M;
% K2 = 1/2*imag(Q_m'*(J_2-J_1)*Q_M);
KK = Q_m'*(J_2)*Q_M;
K1 = real(KK);
K2 = imag(KK);

XX = K1*E_s;
YY = K2*E_s;

% invariance equation %K1*E_s*Upsilon= K2*E_s;
Upsilon = estimation_TLS(XX,YY,K);

% Find the eigenvalues of Psi
Omega = eig(Upsilon);
mu = 2*atan(Omega);
theta_r = asin(-mu/pi);
est_angle = real(rad2deg(theta_r));
end