%% Programmed by ZZF: 2021.12.10
% Input:
% Rxx:          The received sample covariance matrix
% lambda:       Wavelength
% M0:           Total number of antenna elements in Rx array
% K:            Number of signal sources
% L:            The number of subarrays
% L_fb:         The number of subarrays
% m:            The total number of antenna elements in each subarray
% f1:           The handle of the figure
%
% Output:
%

function DOA_run(h_Rxx,lambda,M0,K,L,L_fb,m)

%% Various DOA estimation methods
angle_scan = -90:0.01:90; % Scanning angle range
% [P_Bartlett]=Bartlett_DOA(h_Rxx,M0,lambda,angle_scan);
% [P_LP]=LP_DOA(h_Rxx,M0,lambda,angle_scan);
% [P_ML]=ML_DOA(h_Rxx,M0,lambda,angle_scan);
% [P_capon]=Capon_DOA(h_Rxx,M0,lambda,angle_scan);
[P_MUSIC]=MUSIC_DOA(h_Rxx,M0,lambda,angle_scan,K);
[P_MUSIC_SS]=MUSIC_SS_DOA(h_Rxx,M0,L,lambda,angle_scan,K);
[P_MUSIC_FBSS]=MUSIC_FBSS_DOA(h_Rxx,M0,L_fb,lambda,angle_scan,K);
[est_angle]=ESPRIT_DOA(h_Rxx,m,K);
% [est_angle]=ESPRIT_SS_DOA(h_Rxx,M0,m,K,L);

indd1 = find(P_MUSIC==max(P_MUSIC));
% P_MUSIC_SS(indd1) = 0;
indd2 = find(P_MUSIC_FBSS==max(P_MUSIC_FBSS));
% 
est_ang1 = angle_scan(indd1);
est_ang2 = angle_scan(indd2);


%% Display the result
% plot(angle_scan,10*log10(P_Bartlett/max(P_Bartlett)));

% plot(angle_scan,10*log10(P_capon/max(P_capon)),'Linewidth',1);
% plot(angle_scan,10*log10(P_LP/max(P_LP)),'Linewidth',1);
plot(angle_scan,10*log10(P_MUSIC/max(P_MUSIC)),'Linewidth',1);
grid on
hold on
stem(est_angle,zeros(length(est_angle),1),'or','Linewidth',2);
plot(angle_scan,10*log10(P_MUSIC_SS/max(P_MUSIC_SS)),'Linewidth',1);
plot(angle_scan,10*log10(P_MUSIC_FBSS/max(P_MUSIC_FBSS)),'Linewidth',1);



% for i = 1:K
    msg = num2str(est_angle);
    xline(est_angle,'--r',msg,'LabelOrientation','horizontal',...
        'LabelHorizontalAlignment','center',...
        'LabelVerticalAlignment','top')
    msg = num2str(est_ang1);
    xline(est_ang1,'--r',msg,'LabelOrientation','horizontal',...
        'LabelHorizontalAlignment','center',...
        'LabelVerticalAlignment','middle')
        msg = num2str(est_ang2);
    xline(est_ang2,'--r',msg,'LabelOrientation','horizontal',...
        'LabelHorizontalAlignment','center',...
        'LabelVerticalAlignment','bottom')
% end

% legend('Bartlett','MVDR (Capon)','Linear Prediction','MUSIC','ESPRIT','MUSIC with SS','MUSIC with FBSS')
legend('MUSIC','ESPRIT','MUSIC with SS','MUSIC with FBSS')
% legend('MUSIC','ESPRIT','MUSIC with SS')


title('DOA Estimation');
xlabel('Angle (degree)');
ylabel('Spatila Power Spectrum |P| (dB)');