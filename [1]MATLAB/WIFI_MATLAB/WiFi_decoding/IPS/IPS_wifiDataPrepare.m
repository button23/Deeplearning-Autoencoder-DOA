%% Initialization
% clear;
% clc;
close all;

%% wlan waveform
cfgVHT = wlanVHTConfig('ChannelBandwidth','CBW20','MCS', 0); % 0: BPSK

% Scaling for USRP input
scalingFactor = 2^13;

% Default value is 1024
cfgVHT.APEPLength = 1;

% Create random bit stream
rng(25)
PSDU = randi([0 1],cfgVHT.PSDULength*8,1);

% Waveform generation
wholeWaveform = wlanWaveformGenerator(PSDU,cfgVHT);

% Scaling
wholeWaveform_scaling = scalingFactor * wholeWaveform;
