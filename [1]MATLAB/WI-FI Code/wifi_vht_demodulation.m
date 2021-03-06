%% Demodulate WIFI VHT frame
% Input:
% cfgVHT: Configuration of WIFI VHT frame
% Output:
% chanEst:  Channel state information
% NOTE: The purpose of this code is to create VHT packet on the
% receive end and decodes the received packets received by the USRPs.
%
% Created by ZZF
% Date: March 9, 2022
function [eqSym, chanEst] = wifi_vht_demodulation(cfgVHT,rx)
ind = wlanFieldIndices(cfgVHT);
fs = wlanSampleRate(cfgVHT);       % Sampling rate

% Packet detect and determine coarse packet offset
coarsePktOffset = wlanPacketDetect(rx,cfgVHT.ChannelBandwidth);
if isempty(coarsePktOffset) % If empty no L-STF detected; packet error
    disp('No packet detected');
end

% Extract L-STF and perform coarse frequency offset correction
lstf = rx(coarsePktOffset+(ind.LSTF(1):ind.LSTF(2)),:);
coarseFreqOff = wlanCoarseCFOEstimate(lstf,cfgVHT.ChannelBandwidth);
rx = helperFrequencyOffset(rx,fs,-coarseFreqOff);
% 
% % Extract the non-HT fields and determine fine packet offset
nonhtfields = rx(coarsePktOffset+(ind.LSTF(1):ind.LSIG(2)),:);
finePktOffset = wlanSymbolTimingEstimate(nonhtfields,...
    cfgVHT.ChannelBandwidth);

% Determine final packet offset
pktOffset = coarsePktOffset+finePktOffset;

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

% figure
% plot(20*log10(abs(chanEst)));
% grid on;
% title('Estimated Channel Response');
% xlabel('Subcarrier index');
% ylabel('Power (dB)');

% Extract VHT Data samples from the waveform
vhtdata = rx(pktOffset+(ind.VHTData(1):ind.VHTData(2)),:);

% Estimate the noise power in VHT data field
nVarVHT = vhtNoiseEstimate(vhtdata,chanEstSSPilots,cfgVHT);

% Recover the transmitted PSDU in VHT Data
[~,~,eqSym]  = wlanVHTDataRecover(vhtdata,chanEst,nVarVHT,cfgVHT);
% scatterplot(reshape(eqSym,[],1))
% Determine if any bits are in error, i.e. a packet error
%     biterror = biterr(txPSDU,rxPSDU);
end