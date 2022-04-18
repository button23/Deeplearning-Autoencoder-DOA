% Configuration of WiFi signal and the channel
% DATE: 2022.04.06 by ZZF
%
% Input:
%   cfgVHT:                The configuration of WiFi signal
%   snr:                Signal noise ratio (can be a vector)
%   flag_channel_on:    If add channel 
% Output:
%   cfgVHT:             WiFi VHT signal configuration
%   tgacChannel:        Channel configuration
function [cfgVHT, tgacChannel]=Configurations()
% Create a format configuration object for a SISO VHT transmission
cfgVHT = wlanVHTConfig;
cfgVHT.ChannelBandwidth = 'CBW20'; % Transmitted signal bandwidth
cfgVHT.NumTransmitAntennas = 1;    % Transmit antennas
cfgVHT.NumSpaceTimeStreams = 1;    % Space-time streams
cfgVHT.APEPLength = 0;          % APEP length in bytes
cfgVHT.MCS = 3;                    % Single spatial stream, 64-QAM

%% Channel Impairments
% Parameterize the channel
fs = wlanSampleRate(cfgVHT);       % Sampling rate

tgacChannel = wlanTGacChannel;
tgacChannel.DelayProfile = 'Model-B';
tgacChannel.NumTransmitAntennas = cfgVHT.NumTransmitAntennas;
tgacChannel.NumReceiveAntennas = 1;
tgacChannel.LargeScaleFadingEffect = 'None';
tgacChannel.ChannelBandwidth = cfgVHT.ChannelBandwidth;
tgacChannel.TransmitReceiveDistance = 5;
tgacChannel.SampleRate = fs;
tgacChannel.RandomStream = 'Global stream'; % Global stream  mt19937ar with seed
% tgacChannel.Seed = 10;
end
