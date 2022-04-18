% Compute the Total Least Squares (TLS) estimation
% DATE: 2021.12.10 by ZZF
%
% Input:
%   X:     The sample matrix
%   Y:     The result matrix
%   K:     The number of signal sources

% Output:
%   est_angle:  Estimated angles


function [est_TLS]=estimation_TLS(X,Y,K)

z = [X Y]; % put matrix x and vector y together
[~ , ~, V_z]=svd(z);

est_TLS = -V_z(1:K,K+1:end)*inv(V_z(K+1:end,K+1:end)); % the most important line!!!
% xtyt = -z*V(:,n+1)*V(:,n+1)'; % Find the error of x and y (x_tilda and y_tilda)
%
% % Real Line fitting with the estimated parameter a_tls
% y_tls2 = (XX+xtyt(:,1:n))*a_tls; % (x+x_tilda) multiply by the parameter a_tls
end