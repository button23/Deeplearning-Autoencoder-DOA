% NOTE: The purpose of this code is to create VHT packet on the
% receive end and transmitted by USRP
%
% Created by ZZF
% Date: March 9, 2022
%% Tx USRP Initialization
% release(radio);
radio = comm.SDRuTransmitter(...
    'Platform',             'X300', ...
    'IPAddress',            '192.168.1.24', ...
    'MasterClockRate',      200e6, ...
    'CenterFrequency',      2.4e9, ...
    'Gain',                 20, ...
    'InterpolationFactor',  500, ...
    'ClockSource',          'External', ...
    'PPSSource',            'External', ...
    'ChannelMapping',       [1 2]);

%% Create a format configuration object for a SISO VHT transmission
cfgVHT = wlanVHTConfig;
cfgVHT.ChannelBandwidth = 'CBW40'; % Transmitted signal bandwidth
cfgVHT.NumTransmitAntennas = 1;    % Transmit antennas
cfgVHT.NumSpaceTimeStreams = 1;    % Space-time streams
cfgVHT.APEPLength = 2048;          % APEP length in bytes
cfgVHT.MCS = 5;                    % Single spatial stream, 64-QAM
fs = wlanSampleRate(cfgVHT);       % Sampling rate


% Display the spectrum of the transmitted and received signals. The
% received signal spectrum is affected by the channel
spectrumAnalyzer  = dsp.SpectrumAnalyzer('SampleRate',fs, ...
    'AveragingMethod','Exponential','ForgettingFactor',0.99, ...
    'YLimits',[-30 10],'ShowLegend',true, ...
    'ChannelNames',{'Transmitted waveform'});

% Generate the VHT frame (IEEE 802.11ac)
tx = wifi_vht_waveform_generate(cfgVHT);

%% Transmition starts !
spectrumAnalyzer(tx)

% Inside a while loop, transmit the sine wave using the tx System object.
% Display a message when transmission is complete. Release the radio System object.
while true
    radio(tx);  %data  send_data_scale
end
