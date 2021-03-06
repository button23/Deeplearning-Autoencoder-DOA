% NOTE: The purpose of this code is to receive VHT packet using USRP
% and receive the channel state information

% Created by ZZF
% Date: March 9, 2022
%% Multiple Rx USRPs Initialization
% close all;clc
% release(spectrumScope);
% release(radio);
%     'IPAddress',            '192.168.1.20,192.168.1.21,192.168.1.22,192.168.1.23', ...

radio = comm.SDRuReceiver(...
    'Platform',             'X300', ...
    'IPAddress',            '192.168.1.23', ...
    'MasterClockRate',      200e6, ...
    'CenterFrequency',      2.4e9, ...
    'Gain',                 1, ...
    'DecimationFactor',     10, ...
    'SamplesPerFrame',      4000, ...
    'OutputDataType',       'double', ...
    'ClockSource',          'External', ...
    'PPSSource',            'External', ...
    'ChannelMapping',       [1] );

%% Rx antenna setting
cf = 2.4e9; % 10e6 2.45e9
lambda = physconst('LightSpeed') / cf;
M0 = 8; % Total number of antenna elements
nsnapshot = radio.SamplesPerFrame; % number of snapshots

K = 1; % # of source signal
L = 4; % # of subarrays for FSS
L_fb = 4; % # of subarrays for FBSS   K/2
m = M0-1; % # of antenna elements in ESPRIT subarrays

%% set up time and spectraul scope

frameduration = radio.SamplesPerFrame/(radio.MasterClockRate/radio.DecimationFactor);
timeScope = timescope('TimeSpanSource','Property','TimeSpan',...
    10000*frameduration,'SampleRate',radio.MasterClockRate/500);
% spectrumScope = dsp.SpectrumAnalyzer('SampleRate',200e6/500);
% spectrumScope.ReducePlotRate = true;
%     spectrumScope(NormalizedData);
disp("Reception Started");

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
    'YLimits',[-30 10],'ShowLegend',true, ...
    'ChannelNames',{'Received waveform'});


%% Flag setting
close all
f_normalize = 0;
f_phaseEst = 0;
closeDOAPlot = 1;
f_doa = 0;

countScope = 0;
f1 = figure;
% Receiver start !
% Inside a while loop, receive the sine wave using the rx System object.
% Normalize the signal with respect to amplitude for each receive channel.
% Compute the fast fourier transform (FFT) of each normalized signal100.
% Calculate the phase difference between channel 1 and channel 2, channel 1 and channel 3, and channel 1 and channel 4.

NUM_SAMPLE = 100;
% dataset_wifi = zeros(NUM_SAMPLE,radio.SamplesPerFrame,length(radio.ChannelMapping));
dataset_wifi = zeros(4000*NUM_SAMPLE,1);
save('dataset_wifi.mat','dataset_wifi')
% phaseCompensatedData = zeros(nsnapshot, M0);
% for num = 1 : NUM_SAMPLE
num = 1;
while 1

    [data,len] = step(radio);

    spectrumAnalyzer(data)

    if len == radio.SamplesPerFrame
        dataset_wifi((num-1)*len+1:num*len,:) = data;
    end
    num = num + 1;

    if num == 101
        num = 1;
    end   


    tx = dataset_wifi(200000:300000);
    figure
    plot(real(dataset_wifi(200000:12000+200000)))
    rxx = dataset_wifi(200000:12000+200000);
% 
    rx = wifi_vht_demodulation(cfgVHT,rxx);
    scatterplot(reshape(dataset_wifi,[],1))
%     if len > 0
%         if f_phaseEst
%             [estimatedPhaseOffset]=phase_corr(data,f_normalize);
%             save('estimatedPhaseOffset','estimatedPhaseOffset')
%         end
%         load('estimatedPhaseOffset.mat')
%         % phase offset compensation
%         phaseCompensatedData = data;
%         for i = 1:M0-1
%             phaseCompensatedData(:,i+1) = estimatedPhaseOffset{i}(data(:,i+1));
%         end
    timeScope(real(data));
%         % Enable power measurement.
%         if mod(countScope,100) == 0
%             closeDOAPlot = 0;
%             clf(f1)
%         end
%         %% Doing DOA
%         if ~closeDOAPlot
%             closeDOAPlot = 1;
%             disp('DOA estimation!')
%             %%ind the spacial covariance matrix,Rxx, of the received signal
%             % use the sample average hat{Rxx} to estimate the Rxx
%             phaseCompensatedData_r = phaseCompensatedData.';
%             h_Rxx = phaseCompensatedData_r*phaseCompensatedData_r'/nsnapshot;
%             
%             DOA_run(h_Rxx,lambda,M0,K,L,L_fb,m);
%         end
%     end
%     countScope = countScope + 1;
end
% release(timeScope);
% release(spectrumScope);
% release(radio);
