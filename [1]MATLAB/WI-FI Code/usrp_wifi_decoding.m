%% VHT Configuration
% Create a format configuration object for a SISO VHT transmission
cfgVHT = wlanVHTConfig;
cfgVHT.ChannelBandwidth = 'CBW20'; % Transmitted signal bandwidth
cfgVHT.NumTransmitAntennas = 1;    % Transmit antennas
cfgVHT.NumSpaceTimeStreams = 1;    % Space-time streams
cfgVHT.APEPLength = 1024;          % APEP length in bytes
cfgVHT.MCS = 5;                    % Single spatial stream, 64-QAM
fs = wlanSampleRate(cfgVHT);       % Sampling rate


% Display the spectrum of the transmitted and received signals. The
% received signal spectrum is affected by the channel
spectrumAnalyzer  = dsp.SpectrumAnalyzer('SampleRate',fs, ...
    'AveragingMethod','Exponential','ForgettingFactor',0.99, ...
    'YLimits',[50 100],'ShowLegend',true, ...
    'ChannelNames',{'Received waveform'});


%% 
input = reshape(VarName1,2,length(VarName1)/2);
dataInput = input(2,:) + 1i*input(1,:);
tx = dataInput.';
spectrumAnalyzer(tx)

%%
% rx = wifi_vht_demodulation(cfgVHT,tx);
ind = wlanFieldIndices(cfgVHT);
% fs = wlanSampleRate(cfgVHT);       % Sampling rate
% fs = 20*10^6;
rx = tx(1:4500);
% Packet detect and determine coarse packet offset
coarsePktOffset = wlanPacketDetect(rx,cfgVHT.ChannelBandwidth);
if isempty(coarsePktOffset) % If empty no L-STF detected; packet error
    disp('No packet detected');
end

%%
figure
plot(real(rx))
hold on 
msg = num2str(coarsePktOffset);
xline(coarsePktOffset,'--r',msg,'LabelOrientation','horizontal',...
    'LabelHorizontalAlignment','center',...
    'LabelVerticalAlignment','top')


%%
% Extract L-STF and perform coarse frequency offset correction
lstf = rx(coarsePktOffset+(ind.LSTF(1):ind.LSTF(2)),:);
figure
plot(abs(lstf))
coarseFreqOff = wlanCoarseCFOEstimate(lstf,cfgVHT.ChannelBandwidth);
rx = helperFrequencyOffset(rx,fs,-coarseFreqOff);
figure
plot(real(rx))
%%
close all
figure
plot(real(rx(241:330)))


%%
% 
% % Extract the non-HT fields and determine fine packet offset
nonhtfields = rx(coarsePktOffset+(ind.LSTF(1):ind.LSIG(2)),:);
finePktOffset = wlanSymbolTimingEstimate(nonhtfields,...
    cfgVHT.ChannelBandwidth);

% Determine final packet offset
pktOffset = coarsePktOffset+finePktOffset;
figure
plot(abs(rx(pktOffset + 1:end)))
% Extract L-LTF and perform fine frequency offset correction
lltf = rx(pktOffset+(ind.LLTF(1):ind.LLTF(2)),:);
fineFreqOff = wlanFineCFOEstimate(lltf,cfgVHT.ChannelBandwidth);
rx = helperFrequencyOffset(rx,fs,-fineFreqOff);


%%
close all
for i = 400 : 450
rx_ = rx(i:end);

% Extract VHT-LTF samples from the waveform, demodulate and perform
% channel estimation
vhtltf = rx_((ind.VHTLTF(1):ind.VHTLTF(2)),:);
vhtltfDemod = wlanVHTLTFDemodulate(vhtltf,cfgVHT);

% Get single stream channel estimate
chanEstSSPilots = vhtSingleStreamChannelEstimate(vhtltfDemod,cfgVHT);

% Channel estimate
chanEst = wlanVHTLTFChannelEstimate(vhtltfDemod,cfgVHT);

% figure
% plot(20*log10(abs(chanEst)));
% grid on;
% title('Estimated Channel Response');
% xlabel('Subcarrier index');
% ylabel('Power (dB)');

% Extract VHT Data samples from the waveform
vhtdata = rx_((ind.VHTData(1):ind.VHTData(2)),:);

% Estimate the noise power in VHT data field
nVarVHT = vhtNoiseEstimate(vhtdata,chanEstSSPilots,cfgVHT);

% Recover the transmitted PSDU in VHT Data
[~,~,eqSym]  = wlanVHTDataRecover(vhtdata,chanEst,nVarVHT,cfgVHT);
scatterplot(reshape(eqSym,[],1))
% Determine if any bits are in error, i.e. a packet error
%     biterror = biterr(txPSDU,rxPSDU);
end