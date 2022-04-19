%% Tx antenna setting
close all

nvar = 0.1; % noise power 0.01
method = 'frequency domain';
cf = 5.8e9; % 10e6 2.45e9
lambda = physconst('LightSpeed') / cf;
M0 = 8; % Total number of antenna elements
nsnapshot = 24; % for wifi signals, only 80, 160, ... 800
K = 1; % # of source signal
angle = [-60,-30,20,50,50,70,-30,-45];
coherent_on = 1; % if make the sources coherent to each other
L = K; % # of subarrays for FSS
L_fb = K; % # of subarrays for FBSS   K/2
M = M0 - L + 1; % # of antenna elements in a subarray
m = M0-1; % # of antenna elements in ESPRIT subarrays
msg = ['The SNR is ',num2str(1/nvar), 'dB'];
disp(msg)


%%
% rxArraySize = [4 1];
% 
% rxArray = arrayConfig("Size", rxArraySize, "ElementSpacing", lambda/2);
% 
% % antPosSTA = [repmat(kron(xGridSTA, ones(1, length(yGridSTA))), 1, length(zGridSTA)); ...
% %           repmat(yGridSTA, 1, length(xGridSTA)*length(zGridSTA)); ...
% %           kron(zGridSTA, ones(1, length(yGridSTA)*length(xGridSTA)))];
% antPosSTA = [0;0;0];
% 
% STAs = rxsite("geographic", ...
%     "Antenna", rxArray, ...
%     "AntennaPosition", antPosSTA, ...
%     "AntennaAngle", [0;90], ...
%     "Name", 'a');

% show(STAs)

%% Extract antenna data
loc = 3;
AP1 = squeeze(features(:,:,1,:));
AP2 = squeeze(features(:,:,2,:));
AP3 = squeeze(features(:,:,3,:));
AP4 = squeeze(features(:,:,4,:));

location2_1 = squeeze(AP1(:,:,loc)).'; 
location2_2 = squeeze(AP2(:,:,loc)).'; 
location2_3 = squeeze(AP3(:,:,loc)).'; 
location2_4 = squeeze(AP4(:,:,loc)).'; 
% 
figure
plot(real(location2_1(3,:)))
%% Predict the DOA using MUSIC
% Find the spacial covariance matrix,Rxx, of the received signal
% use the sample average hat{Rxx} to estimate the Rxx
h_Rxx_2_1 = location2_1*location2_1'/24; % 24 samples for one CSI
h_Rxx_2_2 = location2_2*location2_2'/24; % 24 samples for one CSI
h_Rxx_2_3 = location2_3*location2_3'/24; % 24 samples for one CSI
h_Rxx_2_4 = location2_4*location2_4'/24; % 24 samples for one CSI

h_Rxx = double(h_Rxx_2_1);



%%
% doa = zeros(3871,4);
% for i = 1 :3871
%     for j = 1 :4
%         mid = rays{j,i};
% 
%         if isempty(mid)
%             doa(i,j) = 0;
%         else
%             doa(i,j) = mid(1).AngleOfArrival(1);
%         end
%     end
% end
% %% add 180 to make all the angles to be positive
% doa_shift_180 = doa +180;
% doa_shift_180_reshape = reshape(doa_shift_180,3871*4,1);
% min(doa_shift_180_reshape)
% doa_input = sum(doa_shift_180, 2);
% doa_input = int16(doa_input);
% doa_input = doa_input - min(doa_input);
% 
% uniquevalue = unique(doa_input);
% save('doa_input.mat','doa_input')
% %%
% min(doa_input)