% Decoder the WiFi signals and extract the CSI from the WiFi signals
% DATE: 2022.04.06 by ZZF
%
% Input:
%   cfg:            The configuration of WiFi signal
%   tx:             The received WiFi tx signal
%   flag_cfo_est:   The flag for doing CFO estimation
%   plot_on:        Plot the constelation map
%   nsnapshot:      The number of samples used for DOA (ONLY 80,160,...,800 possible)

% Output:
%   rx:             Generated WiFi signal

function [rx,chanEst] = wifi_decoder(cfg,tx,flag_cfo_est,plot_on,nsnapshot)

ind = wlanFieldIndices(cfg);
fs = wlanSampleRate(cfg);       % Sampling rate

% Packet detect and determine coarse packet offset
coarsePktOffset = wlanPacketDetect(tx,cfg.ChannelBandwidth);
if isempty(coarsePktOffset) % If empty no L-STF detected; packet error
    disp('No packet detected');
end

% Extract L-STF and perform coarse frequency offset correction
lstf = tx(coarsePktOffset+(ind.LSTF(1):ind.LSTF(2)),:);

if flag_cfo_est
    coarseFreqOff = wlanCoarseCFOEstimate(lstf,cfg.ChannelBandwidth);
    rx = helperFrequencyOffset(tx,fs,-coarseFreqOff);
end

% Extract the non-HT fields and determine fine packet offset
nonhtfields = rx(coarsePktOffset+(ind.LSTF(1):ind.LSIG(2)),:);
finePktOffset = wlanSymbolTimingEstimate(nonhtfields,...
    cfg.ChannelBandwidth);

% Determine final packet offset
pktOffset = coarsePktOffset+finePktOffset;

% If packet detected outwith the range of expected delays from the
% channel modeling; packet error
if pktOffset>50
    disp('packet detected outwith the range of expected delays')
end

% Extract L-LTF and perform fine frequency offset correction
lltf = rx(pktOffset+(ind.LLTF(1):ind.LLTF(2)),:);

if flag_cfo_est
    fineFreqOff = wlanFineCFOEstimate(lltf,cfg.ChannelBandwidth);
    rx = helperFrequencyOffset(rx,fs,-fineFreqOff);
end

% Extract VHT-LTF samples from the waveform, demodulate and perform
% channel estimation
vhtltf = rx(pktOffset+(ind.VHTLTF(1):ind.VHTLTF(2)),:);
vhtltfDemod = wlanVHTLTFDemodulate(vhtltf,cfg);


% Channel estimate
chanEst = wlanVHTLTFChannelEstimate(vhtltfDemod,cfg);
if plot_on
    figure
    plot(20*log10(abs(chanEst)));
    grid on;
    title('Estimated Channel Response');
    xlabel('Subcarrier index');
    ylabel('Power (dB)');
end
rx = rx(pktOffset+1:pktOffset+nsnapshot,:);

% % Get single stream channel estimate
% chanEstSSPilots = vhtSingleStreamChannelEstimate(vhtltfDemod,cfg);
% 
% Extract VHT Data samples from the waveform
% vhtdata = rx(pktOffset+(ind.VHTData(1):ind.VHTData(2)),:);
% 
% Estimate the noise power in VHT data field
% nVarVHT = vhtNoiseEstimate(vhtdata,chanEstSSPilots,cfg);
% 
% % Recover the transmitted PSDU in VHT Data
% [~,~,eqSym]  = wlanVHTDataRecover(vhtdata,chanEst,nVarVHT,cfg);
% if plot_on
% %     scatterplot(reshape(eqSym,[],1))
% end

% Determine if any bits are in error, i.e. a packet error
%     biterror = biterr(txPSDU,rxPSDU);
end

