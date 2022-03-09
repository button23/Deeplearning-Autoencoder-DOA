%% This code is based on the "Basic WLAN Link Modeling" example code provided by
% Mathworks. The purpose of this code is to create VHT packet on the
% receive end and decodes the received packets received by the USRPs.
% Created by zzf
% Date: March 7, 2022 % edited: March 9, 2022
clear;clc;close all
%% Create a format configuration object for a SISO VHT transmission
cfgVHT = wlanVHTConfig;
cfgVHT.ChannelBandwidth = 'CBW40'; % Transmitted signal bandwidth
cfgVHT.NumTransmitAntennas = 1;    % Transmit antennas
cfgVHT.NumSpaceTimeStreams = 1;    % Space-time streams
cfgVHT.APEPLength = 2048;          % APEP length in bytes
cfgVHT.MCS = 5;                    % Single spatial stream, 64-QAM
fs = wlanSampleRate(cfgVHT);       % Sampling rate

%% Generate the PLCP Preamble and PLCP Header
lstf = wlanLSTF(cfgVHT);
lltf = wlanLLTF(cfgVHT);
lsig = wlanLSIG(cfgVHT);

nonHTfield = [lstf;lltf;lsig]; % Combine the non-HT preamble fields

vhtsiga = wlanVHTSIGA(cfgVHT);
vhtstf = wlanVHTSTF(cfgVHT);
vhtltf = wlanVHTLTF(cfgVHT);
vhtsigb = wlanVHTSIGB(cfgVHT);

preamble = [lstf;lltf;lsig;vhtsiga;vhtstf;vhtltf;vhtsigb];

ind = wlanFieldIndices(cfgVHT);

%% Generate the time-domain VHT data field.
% rng(0) % Initialize the random number generator
% txPSDU = randi([0 1],cfgVHT.PSDULength*8,1); % Generate PSDU data in bits
% Add trailing zeros to allow for channel delay
% txPSDU = [txPSDU; zeros(50,cfgVHT.NumTransmitAntennas)]; %#ok<AGROW>

% data = wlanVHTData(txPSDU,cfgVHT); % txPSDU
% % A VHT waveform is constructed by prepending the non-HT and VHT
% % preamble fields with data
% txWaveform = [preamble;data]; % Transmit VHT PPDU

%% The above code can be replaced using the following code
% txPSDU = [1;0;0;1];
% txWaveform = wlanWaveformGenerator(txPSDU,cfgVHT); % txPSDU

% close all
% figure(1)
% plot(real(txWaveform(1:2000)))
% hold on
% plot(real(txWaveform2(1:2000)))
% legend('without windowing','With windowing')

%% Channel Impairments
% Parameterize the channel
tgacChannel = wlanTGacChannel;
tgacChannel.DelayProfile = 'Model-B';
tgacChannel.NumTransmitAntennas = cfgVHT.NumTransmitAntennas;
tgacChannel.NumReceiveAntennas = 1;
tgacChannel.LargeScaleFadingEffect = 'None';
tgacChannel.ChannelBandwidth = cfgVHT.ChannelBandwidth;
tgacChannel.TransmitReceiveDistance = 5;
tgacChannel.SampleRate = fs;
tgacChannel.RandomStream = 'mt19937ar with seed';
tgacChannel.Seed = 10;

% Display the spectrum of the transmitted and received signals. The
% received signal spectrum is affected by the channel
spectrumAnalyzer  = dsp.SpectrumAnalyzer('SampleRate',fs, ...
    'AveragingMethod','Exponential','ForgettingFactor',0.99, ...
    'YLimits',[-30 10],'ShowLegend',true, ...
    'ChannelNames',{'Transmitted waveform','Received waveform'});

%% Procesing Arrival packets
% 1.The packet is detected.
% 2.Coarse carrier frequency offset is estimated and corrected.
% 3.Fine timing synchronization is established. The L-STF, L-LTF and L-SIG samples are provided
%   for fine timing to allow for packet detection at the start or end of the L-STF.
% 4.Fine carrier frequency offset is estimated and corrected.
% 5.The VHT-LTF is extracted from the synchronized received waveform.
%   The VHT-LTF is OFDM demodulated and channel estimation is performed.
% 6.The VHT Data field is extracted from the synchronized received waveform.
%   The PSDU is recovered using the extracted field and the channel estimate.
% snr = 40:5:50;
snr = 40;
S = numel(snr);
packetErrorRate = zeros(S,1);
% parfor i = 1:S % Use 'parfor' to speed up the simulation
for i = 1:S     % Use 'for' to debug the simulation
    rng(0) % Initialize the random number generator
    
    % Generate a packet waveform
    txPSDU = randi([0 1],cfgVHT.PSDULength*8,1); % PSDULength in bytes
    tx = wlanWaveformGenerator(txPSDU,cfgVHT);
    
    % Add trailing zeros to allow for channel delay
    tx = [tx; zeros(50,cfgVHT.NumTransmitAntennas)]; %#ok<AGROW>
    
    % Pass the waveform through the fading channel model
    reset(tgacChannel); % Reset channel for different realization
    rx = tgacChannel(tx);
    
    % Add noise
    rx = awgn(rx,snr(i));
    
    spectrumAnalyzer([tx rx]);
    
    % Packet detect and determine coarse packet offset
    coarsePktOffset = wlanPacketDetect(rx,cfgVHT.ChannelBandwidth);
    if isempty(coarsePktOffset) % If empty no L-STF detected; packet error
        disp('No packet detected');
        continue; % Go to next loop iteration
    end
    
    % Extract L-STF and perform coarse frequency offset correction
    lstf = rx(coarsePktOffset+(ind.LSTF(1):ind.LSTF(2)),:);
    coarseFreqOff = wlanCoarseCFOEstimate(lstf,cfgVHT.ChannelBandwidth);
    rx = helperFrequencyOffset(rx,fs,-coarseFreqOff);
    
    % Extract the non-HT fields and determine fine packet offset
    nonhtfields = rx(coarsePktOffset+(ind.LSTF(1):ind.LSIG(2)),:);
    finePktOffset = wlanSymbolTimingEstimate(nonhtfields,...
        cfgVHT.ChannelBandwidth);
    
    % Determine final packet offset
    pktOffset = coarsePktOffset+finePktOffset;
    
    % If packet detected outwith the range of expected delays from the
    % channel modeling; packet error
    if pktOffset>50
        disp('packet detected outwith the range of expected delays')
        continue; % Go to next loop iteration
    end
    
    % Extract L-LTF and perform fine frequency offset correction
    lltf = rx(pktOffset+(ind.LLTF(1):ind.LLTF(2)),:);
    fineFreqOff = wlanFineCFOEstimate(lltf,cfgVHT.ChannelBandwidth);
    rx = helperFrequencyOffset(rx,fs,-fineFreqOff);
    
    % Extract VHT-LTF samples from the waveform, demodulate and perform
    % channel estimation
    vhtltf = rx(pktOffset+(ind.VHTLTF(1):ind.VHTLTF(2)),:);
    vhtltfDemod = wlanVHTLTFDemodulate(vhtltf,cfgVHT);
    
    % Get single stream channel estimate
    chanEstSSPilots = vhtSingleStreamChannelEstimate(vhtltfDemod,cfgVHT);
    
    % Channel estimate
    chanEst = wlanVHTLTFChannelEstimate(vhtltfDemod,cfgVHT);
    
    figure
    plot(20*log10(abs(chanEst)));
    grid on;
    title('Estimated Channel Response');
    xlabel('Subcarrier index');
    ylabel('Power (dB)');
    
    % Extract VHT Data samples from the waveform
    vhtdata = rx(pktOffset+(ind.VHTData(1):ind.VHTData(2)),:);
    
    % Estimate the noise power in VHT data field
    nVarVHT = vhtNoiseEstimate(vhtdata,chanEstSSPilots,cfgVHT);
    
    % Recover the transmitted PSDU in VHT Data
    [rxPSDU,~,eqSym]  = wlanVHTDataRecover(vhtdata,chanEst,nVarVHT,cfgVHT);
    scatterplot(reshape(eqSym,[],1))
    
    % Determine if any bits are in error, i.e. a packet error
    %     biterror = biterr(txPSDU,rxPSDU);
end



% constellationDiagram = comm.ConstellationDiagram;
% constellationDiagram.ReferenceConstellation = wlanReferenceSymbols(cfgVHT);
% % Compare received and reference constellation
% constellationDiagram(reshape(eqSym,[],1));
% constellationDiagram.Title = 'Equalized Data Symbols';