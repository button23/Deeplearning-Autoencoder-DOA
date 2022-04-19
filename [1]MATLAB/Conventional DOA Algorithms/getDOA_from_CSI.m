% Compute the DOA using MUSIC, SS MUSIC and FBSS MUSIC
% DATE: 2022.04.11 by ZZF
% 
% Input:
%   rxChan:       The received WiFi signal 

% Output:
%   anglePeak:    The angle of the peak


function [anglePeak] = getDOA_from_CSI(rxChan)
cf = 5.8e9; % 10e6 2.45e9
lambda = physconst('LightSpeed') / cf;
M0 = 8; % Total number of antenna elements
K = 1; % # of source signal
L = 2; % # of subarrays for FSS defualt is K
L_fb = 2; % # of subarrays for FBSS defualt is  K/2
M = M0 - L + 1; % # of antenna elements in a subarray
m = M0-1; % # of antenna elements in ESPRIT subarrays

%%
rxChan = rxChan.';
h_Rxx = rxChan *rxChan'/ length(rxChan);

%%
% nsig = 1;
% [anglePeak,spec,specang] = musicdoa(h_Rxx,nsig);
% % Plot the MUSIC spectrum.
% 
% % plot(specang,10*log10(spec))
% % xlabel('Arrival Angle (deg)')
% % ylabel('Magnitude (dB)')
% % title('MUSIC Spectrum')
% % grid

%% Various DOA estimation methods
angle_scan = -90:0.01:90; % Scanning angle range
% [P_MUSIC]=MUSIC_DOA(h_Rxx,M0,lambda,angle_scan,K);
[P_MUSIC_SS]=MUSIC_SS_DOA(h_Rxx,M0,L,lambda,angle_scan,K);
% [P_MUSIC_FBSS]=MUSIC_FBSS_DOA(h_Rxx,M0,L_fb,lambda,angle_scan,K);
% [est_angle]=ESPRIT_SS_DOA(h_Rxx,M0,m,K,L);
indPeak = find(P_MUSIC_SS == max(P_MUSIC_SS));
anglePeak = angle_scan(indPeak);

%% Display the result
% figure
% plot(angle_scan,10*log10(P_MUSIC/max(P_MUSIC)),'Linewidth',1);
% % plot(angle_scan,10*log10(P_Bartlett/max(P_Bartlett)));
% grid on
% hold on
% % plot(angle_scan,10*log10(P_capon/max(P_capon)),'Linewidth',1);
% % plot(angle_scan,10*log10(P_LP/max(P_LP)),'Linewidth',1);
% % stem(est_angle,zeros(length(est_angle),1),'or','Linewidth',2);
% plot(angle_scan,10*log10(P_MUSIC_SS/max(P_MUSIC_SS)),'Linewidth',1);
% plot(angle_scan,10*log10(P_MUSIC_FBSS/max(P_MUSIC_FBSS)),'Linewidth',1);
% 
% % for i = 1:K
% %     msg = num2str(angle(i));
% %     xline(angle(i),'--r',msg,'LabelOrientation','horizontal',...
% %         'LabelHorizontalAlignment','center',...
% %         'LabelVerticalAlignment','bottom')
% % end
% 
% % legend('Bartlett','MVDR (Capon)','Linear Prediction','MUSIC','ESPRIT','MUSIC with SS','MUSIC with FBSS')
% legend('MUSIC','MUSIC with SS','MUSIC with FBSS')
% 
% title('DOA Estimation');
% xlabel('Angle (degree)');
% ylabel('Spatila Power Spectrum |P| (dB)');
end

