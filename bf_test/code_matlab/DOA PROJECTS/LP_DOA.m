% Compute the Linear Prediction DOA estimate
% DATE: 2021.12.09 by ZZF
% 
% Input:
%   h_Rxx:        The sample covariance matrix 
%   M0: the       Total number of antenna elements in the array
%   lambda:       The wavelength of the signal
%   angle_scan:   The scanning angles during angle search
% Output:
%   P_LP:         Spatial spectrum over scanning angles

function [P_LP]=LP_DOA(h_Rxx,M0,lambda,angle_scan)
mu = -2*pi/lambda*lambda/2*sind(angle_scan);
n = 0 : M0-1;
% all steering vector: M0 (# antenna element) x N (# of scan angle)
sv_all = exp(1i*(n.'-(M0-1)/2)*mu);

P_LP = zeros(length(angle_scan),1);
u = zeros(M0,1);
u(1)=1;
for i = 1:length(angle_scan)
    sv = sv_all(:,i);
    P_LP(i) = abs(u'*inv(h_Rxx)*u/abs(u'*inv(h_Rxx)*sv)^2);
end
% figure
% plot(angle_scan,10*log10(P_Bartlett));
% grid on
% title('Bartlett DOA');
% xlabel('Angle (degree)');
% ylabel('|P(theta)| (dB)');

end