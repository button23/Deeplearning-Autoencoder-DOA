close all

%% Generate a VHT configuration
% Create a VHT configuration
cfgVHT = wlanVHTConfig;
cfgVHT.ChannelBandwidth = 'CBW20';
cfgVHT.NumTransmitAntennas = 1;
cfgVHT.NumSpaceTimeStreams = 1;
cfgVHT.MCS = 4; %4:16QAM  8:256QAM          % 16-QAM and 3/4 coding rate
cfgVHT.APEPLength = 1024; % Bytes

% Create a random PSDU
s = rng(10); % Seed the random number generator
psdu = randi([0 1],cfgVHT.PSDULength*8,1,'int8');

% Generate a VHT packet
tx = wlanWaveformGenerator(psdu,cfgVHT);

%% Model Impairments
%sampling frequency offset between the transmitter and receiver is modeled
%by resampling the transmitted waveform. The resample function can be used
%to model a limited range of large sample frequency offsets. In this example
%a -100 parts per million (PPM) sample rate offset is modeled.

% Model sample rate offset
p = 1e4;   % Interpolation factor
q = 1e4-1; % Decimation factor
sro = (1-p/q)*1e6;
disp('Impairments:');
disp(['  Sample rate offset (SRO): ' num2str(sro,'%.1f') ' PPM']);

% Resample the waveform, appending zeros to allow for filter delay
N = 100; % Size of filter used for resampling
rx = resample([tx; zeros(N,cfgVHT.NumTransmitAntennas)],p,q,N);
% rx = tx;
%% Add Residual Frequency Offset
fc = 5.25e9;         % Carrier frequency, Hertz
cfo = (sro*1e-6)*fc; % Carrier frequency offset, Hertz
disp(['  Carrier frequency offset (CFO): ' num2str(cfo,'%.1f') ' Hz']);

fs = wlanSampleRate(cfgVHT);           % Baseband sample rate
rx = helperFrequencyOffset(rx,fs,cfo); % Add frequency offset

%% Add noise with 30dBW variance
awgnChannel = comm.AWGNChannel('NoiseMethod','Variance','Variance',10^(-30/10));
rx = awgnChannel(rx);

% rx = synchFrame;

%% Front End Synchronization and Receiver Processing
% Generate field indices
ind = wlanFieldIndices(cfgVHT);

% Packet detection
tOff = wlanPacketDetect(rx,cfgVHT.ChannelBandwidth);

%% Coarse frequency offset correction
lstf = rx(tOff+(ind.LSTF(1):ind.LSTF(2)),:);
coarseCFOEst = wlanCoarseCFOEstimate(lstf,cfgVHT.ChannelBandwidth);
% rx = helperFrequencyOffset(rx,fs,-coarseCFOEst);

%% Self-made: Coarse CFO Estimation
% Coarse CFO Estimation
% Using the last five identical parts to do coarse frequency offset
% estimation
fisrtPartStf = lstf(81:144);
secondPartStf = lstf(97:160);
% Caluculate the phase shift between every sample (unit:radian)
phase_STF = 1/16 * angle(fisrtPartStf'*secondPartStf);
% the calculated phase is in the unit of radian,change it to frequency.
Ts = 1/(2*10^7);
CFO_STF = int32(phase_STF/(2*pi*Ts));  % CFO

% Coarse CFO Correlation
% Compensate only the LTF, and the index is set to be zero instead of
% the initial index in STF to symplify the process
% because LTF would do later calibration any way
mmC = 0 : 127;
kkc = 0:length(rx)-1;   % !!! Changing the index(e.g. starts from a ramdom number) does not have a big effect on the result. I guess the reason is that the CFO is so small that just using pilot can fix it. 
% Be caucious, this line can not be run multiple times, because the rxSig
% would be overwritten
% rxSig(ind_LTF(1)+32:ind_LTF(2)) = rxSig(ind_LTF(1)+32:ind_LTF(2)) .* exp(-1i * mmC' * phase_STF);
rx = rx.* exp(-1i * kkc' * phase_STF);  % Confirmed correct compared to the result !!!!
rx2 = rx.* exp(-1i * (kkc)' * (phase_STF+0.11));  % Confirmed correct compared to the result !!!!
rx3 = rx.* exp(-1i * (kkc+2000)' * (phase_STF));  % Confirmed correct compared to the result !!!!

% rx4 = rx.* exp(-1i * mmC' * (phase_STF));  % Confirmed correct compared to the result !!!!

% BIG QUESTION!
% If I change the scale of the index, it won't hurt the
% performance,however; if I touch on phase_STF, say change add 0.1 to it,
% the performance would be severely damaged. What would be the reason ????
% Attempted answer:
% I guess it is due to the linearity of the index. Even though the start of
% the index may be changed to a big number. However, since it increments by
% 1. ......


%% Symbol timing synchronization (It applies after coarse carrier frequency compensation)
nonhtPreamble = rx(tOff+(ind.LSTF(1):ind.LSIG(2)),:); % Extract LSTF again after this part has been CFO compensated. 
symOff = wlanSymbolTimingEstimate(nonhtPreamble,cfgVHT.ChannelBandwidth);
tOff = tOff+symOff;
% tOff = -14;

%% Fine CFO Estimation
% Fine frequency offset correction
lltf = rx(tOff+(ind.LLTF(1):ind.LLTF(2)),:);
fineCFOEst = wlanFineCFOEstimate(lltf,cfgVHT.ChannelBandwidth);
% rx = helperFrequencyOffset(rx,fs,-fineCFOEst);

%% Self-made code? Fine CFO Estimation
first_LTF = lltf(33:96);
second_LTF = lltf(97:160);

phase_LTF = 1/64 * angle(first_LTF'*second_LTF);
CFO_LTF = phase_LTF/(2 * pi * Ts);
% rx = helperFrequencyOffset(rx,fs,-CFO_LTF);

% % Total carrier frequency offset phase combined with coarse CFO estimation and
% % Fine CFO estimation
total_phase = phase_LTF + phase_STF;
% total_cfo = CFO_LTF + CFO_STF;

% Fine CFO Correlation
mmF = 128:length(rx)-1-128;
% mmF = 128:ind_Data(2);
mmF = double(mmF');

% If do Fine Compensation at the beginning, the error rate is less and the
% constallation map is  neater.
% rx(ind.LSIG(1):end) = rx(ind.LSIG(1):end) .* exp(-j * mmF * phase_LTF);
rx = rx.* exp(-1i * kkc' * phase_LTF);

%% Noise power estimation
lltf = rx(tOff+(ind.LLTF(1):ind.LLTF(2)),:);
lltfDemod = wlanLLTFDemodulate(lltf,cfgVHT);
noiseEst = helperNoiseEstimate(lltfDemod,cfgVHT.ChannelBandwidth,cfgVHT.NumSpaceTimeStreams);
% noiseEst = 0;

%% Channel estimation
vhtltf = rx(tOff+(ind.VHTLTF(1):ind.VHTLTF(2)),:);
vhtltfDemod = wlanVHTLTFDemodulate(vhtltf,cfgVHT);
% chanEst = wlanVHTLTFChannelEstimate(vhtltfDemod,cfgVHT);

%% Self made: Channel Estimation
cp_remove_vhtltf = vhtltf(17:end);
dsp_test_vhtltf = fft(cp_remove_vhtltf);
dsp_test_vhtltf_fftshift = double(fftshift(dsp_test_vhtltf,1))*1/64 *sqrt(56);
vht_llf_af_pilots = [dsp_test_vhtltf_fftshift(5:32);dsp_test_vhtltf_fftshift(34:61)];

seqLong_VHT = [1;1;1;1;-1;-1;1;1;-1;1;-1;1;1;1;1;1;1;
    -1;-1;1;1;-1;1;-1;1;1;1;1;1;-1;-1;1;
    1;-1;1;-1;1;-1;-1;-1;-1;-1;1;1;-1;-1;1;-1;
    1;-1;1;1;1;1;-1;-1];
H_VHT = vht_llf_af_pilots .* seqLong_VHT;
chanEst = H_VHT;

%% Recovery Without Sample Rate Offset or Residual CFO Tracking
disp('Front-end impairment correction:');
frontEndCFOEst = coarseCFOEst+fineCFOEst;
self_EstimatedCFO = CFO_STF + CFO_LTF;
disp(['  Estimated CFO: ' num2str(frontEndCFOEst,'%.1f') ' Hz']);
disp(['  Self Estimated CFO: ' num2str(self_EstimatedCFO,'%.1f') ' Hz']);

residualCFO = cfo-frontEndCFOEst;
self_residualCFO = cfo-self_EstimatedCFO;
disp(['  Residual CFO after initial correction: ' num2str(residualCFO,'%.1f') ' Hz']);
disp(['  Self Residual CFO after initial correction: ' num2str(self_residualCFO,'%.1f') ' Hz']);

%% Recover the VHT data field with optional pilot tracking
%The helper function trackingVHTDataRecover recovers the VHT data field with 
%optional pilot tracking to correct for timing and phase errors due to SRO 
%and CFO. Pilot tracking is controlled using the helper object trackingRecoveryConfig.
% Recovery configuration with pilot tracking disabled
cfgRec = trackingRecoveryConfig;
cfgRec.PilotTracking = 'None';

% Extract data field with Ne additional samples to allow for negative SRO
maxDuration = 5.484e-3; % Maximum packet duration in seconds
maxSRO = 120;           % PPM
Ne = ceil(fs*maxDuration*maxSRO*1e-6); % Number of extra samples
dataInd = tOff+(ind.VHTData(1):ind.VHTData(2)+Ne);
dataInd = dataInd(dataInd<=length(rx)); % Only use indices within waveform
data = rx(dataInd,:);

% Perform demodulation and decoding
[rxPSDUNoTrack,~,eqSymNoTrack] = trackingVHTDataRecover(data,chanEst,noiseEst,cfgVHT,cfgRec);

%% Equalized Constellation
%The equalized constellation is plotted which shows a rotation of all 
%constellation points caused by residual CFO, and a spreading of constellation points due to SRO. 
ConstNoTrack = comm.ConstellationDiagram;
ConstNoTrack.Title = 'Equalized symbols with no pilot tracking';
% ConstNoTrack.ReferenceConstellation = wlanReferenceSymbols(cfgVHT);
ConstNoTrack(eqSymNoTrack(:));

[~,berNoTrack] = biterr(rxPSDUNoTrack,psdu);
disp('Bit error rate:');
disp(['  Without tracking: ' num2str(berNoTrack)]);


%% Recovery With Sample Rate Offset Tracking and Residual CFO Tracking
%Now the data field is recovered with joint timing and phase pilot
%tracking to correct for SRO and residual CFO.

% Recovery configuration with pilot tracking enabled
cfgRec = trackingRecoveryConfig;
cfgRec.PilotTracking = 'Joint'; % Joint timing and phase tracking
cfgRec.PilotTrackingWindow = 9; % Averaging window in OFDM symbols

% Perform demodulation and decoding
[rxPSDU,~,eqSymTrack,cpe,peg] = trackingVHTDataRecover(data,chanEst,noiseEst,cfgVHT,cfgRec);
ConstTrack = comm.ConstellationDiagram;
ConstTrack.Title = 'Equalized symbols with joint pilot tracking';
% ConstTrack.ReferenceConstellation = wlanReferenceSymbols(cfgVHT);
ConstTrack(compensatedPhase(:));

%% Self-made Phase Tracking
residual_offset = zeros(53,1);
ppangle = zeros(3,1);
averPPangle = zeros(53,1);
referenGen = zeros(4,53);
interr_ = zeros(4,53);
pilotFunctt = zeros(4,53);
% zeropadding_remove_origin = eqSymNoTrack;
H_pilot = [H_VHT(8);H_VHT(22);H_VHT(35);H_VHT(49)];

%Data Part
data_part = data(1:53*80);
data_part_re = reshape(data_part,80,53);
data_part_re_de = data_part_re(17:end,:);

[cfgOFDM, dataInd, pilotInd] = wlan.internal.wlanGetOFDMConfig(cfgVHT.ChannelBandwidth, cfgVHT.GuardInterval, 'VHT', 1);
symOffset = 0.75;

postCPRemoval = data_part_re_de * sqrt(56) * 1/64;
postFFT = fft(postCPRemoval, [], 1);
postShift = fftshift(postFFT,1);

zeropadding_remove_origin = [postShift(5:32,:);postShift(34:61,:)];
pilot_data_origin = [zeropadding_remove_origin(8,:);zeropadding_remove_origin(22,:);zeropadding_remove_origin(35,:);zeropadding_remove_origin(49,:)];
% pilot_data_origin = ofdmDemodPilots;

% polarity sequence
% Used for pilot bit (served as sign)
pv=[1;1;1;1; -1;-1;-1;1; -1;-1;-1;-1; 1;1;-1;1; -1;-1;1;1; -1;1;1;-1; 1;1;1;1; 1;1;-1;1;
    1;1;-1;1; 1;-1;-1;1; 1;1;-1;1; -1;-1;-1;1; -1;1;-1;-1; 1;-1;-1;1; 1;1;1;1; -1;-1;1;1;
    -1;-1;1;-1; 1;-1;1;1; -1;-1;-1;1; 1;-1;-1;-1; -1;1;-1;-1; 1;-1;1;1; 1;1;-1;1; -1;1;-1;1;
    -1;-1;-1;-1; -1;1;-1;1; 1;-1;1;-1; 1;1;1;-1; -1;1;-1;-1; -1;1;1;1; -1;-1;-1;-1; -1;-1;-1];

Data_pilot = zeros(64,53);
pil=zeros(4,4);
% The following four sequences are repeatedly used for pilot bit of OFDM
% symbol (the sequence should be multiplied by the corresponding polarity bit)
pil(:,1) = [1,1,1,-1];
pil(:,2) = [1,1,-1,1];
pil(:,3) = [1,-1,1,1];
pil(:,4) = [-1,1,1,1];
for reff = 1:53
    referenGen(:,reff) = pv(reff+4)*pil(:,mod((reff-1),4)+1);
end
for refff = 1:53
    interr = conj(pilot_data_origin(:,refff)).* H_pilot; % 
    interr_(:,refff) = interr.*referenGen(:,refff);
end
% interr_ = perr;
% Averaging
averPilot = zeros(4,53);
omegaa = zeros(53,1);
deltaa = zeros(53,1);
kkpp = [-21,-7,7,21;ones(1,4)]';

for fffsw = 1:53
    % if doing the following averaging, the result becomes worse.
%     if fffsw < 9
%         interer = sum(interr_(:,1:fffsw),2);
%     else
%         interer = sum(interr_(:,(fffsw-8:fffsw)),2);
%     end

    interer = interr_(:,fffsw);
    if fffsw >1
        averPilot(:,fffsw) = interer * exp(-1i*omegaa(fffsw-1));
    else
        averPilot(:,fffsw) = interer;
    end
    jjjs =  lscov(kkpp,angle(averPilot(:,fffsw)),abs(averPilot(:,fffsw)));
    deltaa(fffsw) = jjjs(1);
    if fffsw > 1
       omegaa(fffsw) = jjjs(2)+omegaa(fffsw-1);
   else
       omegaa(fffsw) = jjjs(2);
   end
end

% Experiment (using the calculated phase by Matlab)
phase11 = zeros(52,53);

indexData = -28:28;
indexData(8) = [];
indexData(22) = [];
indexData(29) = [];
indexData(36) = [];
indexData(50) = [];

% phase calculation
for symbbo = 1:53
    phase11(:,symbbo) = exp(1i*deltaa(symbbo)*indexData') * exp(1i*omegaa(symbbo));
end

compensatedPhase = eqSymNoTrack .* phase11;
% compensatedPhase = eqSymNoTrack;

ConstTrack_ = comm.ConstellationDiagram;
ConstTrack_.Title = 'Equalized symbols with joint pilot tracking(zzf)';
% ConstTrack.ReferenceConstellation = wlanReferenceSymbols(cfgVHT);
ConstTrack_(compensatedPhase(:));


% %% Equalized constellation 
% % Shows a clear 16-QAM constellation with no spreading or rotation. There are no bit errors.
% ConstTrack = comm.ConstellationDiagram;
% ConstTrack.Title = 'Equalized symbols with joint pilot tracking';
% % ConstTrack.ReferenceConstellation = wlanReferenceSymbols(cfgVHT);
% ConstTrack(eqSymTrack(:));
% 
% [~,berTrack] = biterr(rxPSDU,psdu);
% disp(['  With tracking: ' num2str(berTrack)]);
% 
% 
% %% The function trackingVHTDataRecover returns measurements from which the residual CFO, and SRO can be estimated
% %cpe - The common phase error (radians) per symbol
% %peg - The phase error gradient (radians per subcarrier) per symbol
% 
% [residualCFOEst,sroEst] = trackingPlotSROCFOEstimates(cpe,peg,cfgVHT);
% 
% % Display estimated SRO, residual CFO and total CFO
% fprintf('Tracked impairments:\n');
% fprintf('  Estimated residual CFO: %3.1f Hz (%.1f Hz error)\n', ...
%     residualCFOEst,residualCFOEst-residualCFO);
% % fprintf('  Estimated residual CFO: %3.1f Hz (%.1f Hz error)\n', ...
% %     residualCFOEst,residualCFOEst-self_residualCFO);
% fprintf('  Estimated SRO: %3.1f PPM (%.1f PPM error)\n',sroEst,sroEst-sro);
% cfoEst = frontEndCFOEst+residualCFOEst; % Initial + tracked CFO estimate
% % cfoEst = self_EstimatedCFO+residualCFOEst; % Initial + tracked CFO estimate
% fprintf('Estimated CFO (initial + tracked): %.1f Hz (%.1f Hz error)\n',cfoEst,cfoEst-cfo);
% 
% rng(s); % Restore the state of the random number generator

