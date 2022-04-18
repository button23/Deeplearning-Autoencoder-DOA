%% Three beamforming algorithms are illustrated. 
% The phase-shift beamformer
% The minimum variance distortionless response (MVDR) beamformer
% The linearly constrained minumum variance (LCMV) beamformer
clear;close all;clc
% Define the incoming signal
t = 0:0.001:0.3;
s = zeros(size(t))';
s(201:205) = 1;
figure
plot(t,s)
title('Pulse');
xlabel('Time (s)');
ylabel('Amplitude (V)');

%% Tx antenna setting
cf = 2.4e9; % 10e6 2.45e9
lambda = physconst('LightSpeed') / cf;
ula = phased.ULA('NumElements',10,'ElementSpacing',lambda/2);
angle = [45;0];

%% Received signal at each array element
rx_ori = collectPlaneWave(ula,s,angle,cf);
% Add thermal noise (complex Gaussian distribution). \
% The power of the noise is 0.5 watt. Here the power of the signal is 1 watt;
% therefore, the SNR is 1/0.5 = 2, which is 3dB.
rs = RandStream.create('mt19937ar','Seed',2008);

nvar = 0.5;
noise = sqrt(nvar/2)*(randn(rs,size(rx_ori)) + 1i*randn(rs,size(rx_ori)));
% from the cov, we can see that the cov(x,x) is 0.5, the cov(x,y) is nearly
% 0. It will be 0 if the size is infinity.
cov = noise'*noise/length(rx_ori); % no use, just to verify.
rx_noise = rx_ori + noise;
% Plot the received signal at first two antenna elements
figure
subplot(211); 
plot(t,abs(rx_noise(:,1)));axis tight;
title('Pulse at Antenna 1');xlabel('Time (s)');ylabel('Magnitude (V)');
subplot(212);
plot(t,abs(rx_noise(:,2)));axis tight;
title('Pulse at Antenna 2');xlabel('Time (s)');ylabel('Magnitude (V)');

%% Phase shift beamformer
psbeamformer = phased.PhaseShiftBeamformer('SensorArray',ula,...
    'OperatingFrequency',cf,'Direction',angle,'WeightsOutputPort',true);
[yPS,wPS]= psbeamformer(rx_noise);
% Plot the output
figure
plot(t,abs(yPS)); axis tight;
title('Output of Phase Shift Beamformer');
xlabel('Time (s)');ylabel('Magnitude (V)');

figure 
% Plot array response with phase-shft weighting for each antenna element
pattern(ula,cf,-180:180,0,'Weights',wPS,'Type','powerdb',...
    'PropagationSpeed',physconst('LightSpeed'),'Normalize',false,...
    'CoordinateSystem','rectangular')
axis([-90 90 -60 0])

%% Interference Signals Modeling
nSamp = length(t);
% 10 means the power is 100, muach stronger than both 
% the tx signal and the noise
s1 = 10*randn(rs,nSamp,1); 
s2 = 10*randn(rs,nSamp,1);
s3 = 0.1*randn(rs,nSamp,1);
s4 = 0.1*randn(rs,nSamp,1);

angle_inter = [70 10 30 25;0 0 0 0];
% Interence signal comming from 30 and 50 degrees in azimuth
interence = collectPlaneWave(ula,[s1 s2 s3 s4],angle_inter,cf);
% Assume the SNR is 50dB so the effect of the noise can be minimized for us
% to observe the effect of the interference signals
noisePwr = 0.00001; %50dB SNR
noise_sml = sqrt(noisePwr/2)*(randn(rs,size(rx_ori)) + 1i*randn(rs,size(rx_ori)));
rx_int = interence + noise_sml;
rx_signal = rx_ori + rx_int;    

%% Phase shift beamformer (failed in case of interference)
yConBF = psbeamformer(rx_signal);
% From the figure, we can see that, because the interference signals are 
% much stronger than the target signal, we cannot extract the signal content.
figure
plot(t,abs(yConBF)); axis tight;
title('Output of Phase Shift Beamformer With Presence of Interference');
xlabel('Time (s)');ylabel('Magnitude (V)');

%% MVDR beamformer (also called the Capon beamformer)
% The MVDR beamformer preserves the signal arriving along a desired 
% direction, while trying to suppress signals coming from other directions. 
% 'TrainingInputPort' means we have access to the interence data
mvdrbeamformer = phased.MVDRBeamformer('SensorArray',ula,...
    'OperatingFrequency',cf,'Direction',angle,...
    'WeightsOutputPort',true,'TrainingInputPort',true);
[yMVDR,wMVDR] = mvdrbeamformer(rx_signal,rx_int);

% Looking at the response pattern of the beamformer, we see two deep nulls 
% along the interference directions, (30 and 50 degrees). 
% The beamformer also has a gain of 0 dB along the target direction 
% of 45 degrees. Thus, the MVDR beamformer preserves the target signal 
% and suppresses the interference signals.  PhaseShift pattern does not 
% null out the interference at al
figure
plot(t,abs(yMVDR)); axis tight;
title('Output of MVDR Beamformer With Presence of Interference');
xlabel('Time (s)');ylabel('Magnitude (V)');

figure
% Plot array response with phase-shft weighting for each antenna element
pattern(ula,cf,-180:180,0,'Weights',wMVDR,'Type','powerdb',...
    'PropagationSpeed',physconst('LightSpeed'),'Normalize',false,...
    'CoordinateSystem','rectangular')
axis([-90 90 -80 20])

hold on;
pattern(ula,cf,-180:180,0,'Weights',wPS,'Type','powerdb',...
    'PropagationSpeed',physconst('LightSpeed'),'Normalize',false,...
    'CoordinateSystem','rectangular')
hold off;
legend('MVDR','Phase Shift')

%% Self-nulling issue in MVDR
% MVDR requires the knowledge of the interence.if the target signal is 
% received along a direction slightly different from the desired one, 
% the MVDR beamformer suppresses it. This occurs because the MVDR beamformer
% treats all the signals, except the one along the desired direction, as 
%undesired interferences.

% To illustrate this self nulling effect, we define an MVDR beamformer 
% and set the TrainingInputPort property to false.

mvdrbeamformer_selfnull = phased.MVDRBeamformer('SensorArray',ula,...
    'Direction',angle,'OperatingFrequency',cf,...
    'WeightsOutputPort',true,'TrainingInputPort',false);

% Set the mismatch in signal direction
expDirection = [43;0]; % The actual direction is 45 degree in azimuth 
mvdrbeamformer_selfnull.Direction = expDirection;
% Apply MVDR beamformer to the received signal
[ySn, wSn] = mvdrbeamformer_selfnull(rx_signal);

% The receiver cannot differentiate the target signal and the interference.
figure
plot(t,abs(ySn)); axis tight;
title('Output of MVDR Beamformer With Signal Direction Mismatch');
xlabel('Time (s)');ylabel('Magnitude (V)');

% Beamformer response pattern
%  the MVDR beamformer tries to suppress the signal arriving along 
% 45 degrees because it is treated like an interference signal.
figure
pattern(ula,cf,-180:180,0,'Weights',wSn,'Type','powerdb',...
    'PropagationSpeed',physconst('LightSpeed'),'Normalize',false,...
    'CoordinateSystem','rectangular');
axis([-90 90 -40 25]);

%% LCMV Beamformer
% LCMV prevents signal self-nulling, it allows up put multiple constrains
% along the target direction (steering vector). 
lcmvbeamformer = phased.LCMVBeamformer('WeightsOutputPort',true);

% The weights are given by the steering vector that steer the array toward
% the directions according to the constrants we want to have on LCMV BF
steeringvec = phased.SteeringVector('SensorArray',ula);
%41 and 45 are +/2 2 degrees of the expected direction (43 degree)
stv = steeringvec(cf,[43 41 45]); 
lcmvbeamformer.Constraint = stv;
lcmvbeamformer.DesiredResponse = [1;1;1];

% Apply LCMV BF to the received signal (including interference)
[yLCMV,wLCMV] = lcmvbeamformer(rx_signal);

% The LCMV response pattern shows that the beamformer puts the constraints 
% along the specified directions, while nulling the interference 
% signals along 30 and 50 degrees.
figure
plot(t,abs(yLCMV)); axis tight;
title('Output of LCMV Beamformer With Signal Direction Mismatch');
xlabel('Time (s)');ylabel('Magnitude (V)');

% Array response pattern between 0 and 90 degrees in azimuth 
figure
pattern(ula,cf,-180:180,0,'Weights',wLCMV,'Type','powerdb',...
    'PropagationSpeed',physconst('LightSpeed'),'Normalize',false,...
    'CoordinateSystem','rectangular');
axis([0 90 -40 35]);
hold on;  % compare to MVDR
pattern(ula,cf,-180:180,0,'Weights',wSn,...
    'PropagationSpeed',physconst('LightSpeed'),'Normalize',false,...
    'Type','powerdb','CoordinateSystem','rectangular');
hold off;
legend('LCMV','MVDR');

%% 2D Array Beamforming 
% Define a URA 
colSp = 0.5*lambda;
rowSp = 0.5*lambda;
ura = phased.URA('Size',[4 4],'ElementSpacing',[rowSp colSp]);

rx_ura = collectPlaneWave(ura,s,angle,cf);
% AWGN
noise = sqrt(nvar/2)*(randn(rs,size(rx_ura)) + 1i*randn(rs,size(rx_ura)));

%% Unlike a ULA, which can only differentiate the angles in azimuth 
%direction, a URA can also differentiate angles in elevation direction. 
%Therefore, we specify two interference signals arriving along the 
% directions [30;10] and [50;-5] degrees, respectively.
ang_inter_ura = [30 50;10 -5];
s1 = 10*randn(rs,nSamp,1);
s2 = 10*randn(rs,nSamp,1);
interference_ura = collectPlaneWave(ura,[s1 s2],ang_inter_ura,cf);
inter_ura = noise + interference_ura;
rx_ura_int = rx_ura + inter_ura;

%% MVDR Beamformer 
mvdrbeamformer_ura = phased.MVDRBeamformer('SensorArray',ura,...
    'OperatingFrequency',cf,'Direction',angle,'WeightsOutputPort',true,...
    'TrainingInputPort',true);
[yMVDR_ura,wMVDR_ura] = mvdrbeamformer_ura(rx_ura_int,inter_ura);

figure
plot(t,abs(yMVDR_ura));axis tight;
title('Output of MVDR Beamformer for URA');
xlabel('Time (s)');ylabel('Magnitude (V)');

% Plot the array response pattern 
subplot(2,1,1);
pattern(ura,cf,-180:180,-5,'Weights',wMVDR_ura,'Type','powerdb',...
    'PropagationSpeed',physconst('LightSpeed'),'Normalize',false,...
    'CoordinateSystem','rectangular');
title('Response Pattern at -5 Degrees Elevation');
axis([-90 90 -60 10]);
subplot(2,1,2);
pattern(ura,cf,-180:180,10,'Weights',wMVDR_ura,'Type','powerdb',...
    'PropagationSpeed',physconst('LightSpeed'),'Normalize',false,...
    'CoordinateSystem','rectangular');
title('Response Pattern at 10 Degrees Elevation');
axis([-90 90 -60 10]);