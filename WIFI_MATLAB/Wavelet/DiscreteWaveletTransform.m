%% Animation version (the window keeps moving)
clear all; close all; clc
L = 10;
n= 2048;
t2 = linspace(0, L, n+1);  t=t2(1:n);
% Target function that we want to analyze
S = (3*sin(2*t) +0.5*tanh(0.5*(t-3))+0.2*exp(-(t-4).^2) +1.5*sin(5*t)+4*cos(3*(t-6).^2))/10 +(t/20).^3;
St = fft(S);
% k = (2*pi/L)*[0:n/2-1 -n/2:-1];   ks = fftshift(k);
ks = (2*pi/L)*(-n/2 : n/2-1);
% the increment is the step with which you move your
% window rightwards and  should be decided by how wide your window is.
tslide = 0 : 0.3 : 10;
sgt_mat = zeros(length(tslide), n);

for j=1:length(tslide)
    % Gabor (Wavelet function)
    g = exp(-(t-tslide(j)).^2); % original window
%     g = exp(-0.05*(t-tslide(j)).^2); % large window
%     g = exp(-5*(t-tslide(j)).^2); % short window
    % Filtering out the outside part of th window
    Sg = g .*S;
    Sgt = fft(Sg);
    sgt_mat(j, :) = abs(fftshift(Sgt));
%     figure(1)
%     subplot(3,1,1),  plot(t, S, 'k', t, g, 'r')
%     subplot(3,1,2),  plot(t, Sg, 'k')
%     subplot(3,1,3),  plot(ks, abs(fftshift(Sgt))/max(abs(fftshift(Sgt))))
%     axis([-50 50 0 1])
%     drawnow %
%     pause(0.3)
end

figure(2)
n = 34
pcolor(tslide(1:n), ks, sgt_mat(1:n,:).'), shading interp
set(gca, 'Ylim', [-50 50])
colormap(hot)
