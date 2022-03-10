%% Generate WIFI VHT frame with random data load
% Input: 
% cfgVHT: Configuration of WIFI VHT frame
% Output:
% tx:     Generated VHT frame given by the input configuration 
% NOTE: The purpose of this code is to create VHT packet on the
% receive end and decodes the received packets received by the USRPs.
%
% Created by ZZF
% Date: March 9, 2022

function tx = wifi_vht_waveform_generate(cfgVHT)

rng(0) % Initialize the random number generator

% Generate a packet waveform
txPSDU = randi([0 1],cfgVHT.PSDULength*8,1); % PSDULength in bytes
tx = wlanWaveformGenerator(txPSDU,cfgVHT);

% Add trailing zeros to allow for channel delay
% tx = [tx; zeros(50,cfgVHT.NumTransmitAntennas)]; %#ok<AGROW>
end