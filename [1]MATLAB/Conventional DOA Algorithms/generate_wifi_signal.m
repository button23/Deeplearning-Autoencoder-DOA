% Generate the WiFi signals with the given SNR
% DATE: 2022.04.06 by ZZF
%
% Input:
%   cfg:                The configuration of WiFi signal
%   cfgChannel:         The configuration of Channel
%   snr:                Signal noise ratio (can be a vector)
%   flag_channel_on:    If add channel 
% Output:
%   rx:             Generated WiFi signal


function [rx]=generate_wifi_signal(cfg,tgacChannel,snr,flag_channel_on)

% Display the spectrum of the transmitted and received signals. The
% received signal spectrum is affected by the channel
% spectrumAnalyzer  = dsp.SpectrumAnalyzer('SampleRate',fs, ...
%     'AveragingMethod','Exponential','ForgettingFactor',0.99, ...
%     'YLimits',[-30 10],'ShowLegend',true, ...
%     'ChannelNames',{'Transmitted waveform','Received waveform'});

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

S = numel(snr);
% parfor i = 1:S % Use 'parfor' to speed up the simulation
for i = 1:S     % Use 'for' to debug the simulation
        % Generate a packet waveform
        txPSDU = randi([0 1],cfg.PSDULength*8,1); % PSDULength in bytes
        tx = wlanWaveformGenerator(txPSDU,cfg);

        % Add trailing zeros to allow for channel delay
        tx = [tx; zeros(50,cfg.NumTransmitAntennas)]; %#ok<AGROW>

        if flag_channel_on
            % Pass the waveform through the fading channel model
            reset(tgacChannel); % Reset channel for different realization
            rx = tgacChannel(tx);
        else
            rx = tx; % No channel
        end
        % Add noise
%         rx = awgn(rx,snr(i));
%     spectrumAnalyzer([tx rx]);
end