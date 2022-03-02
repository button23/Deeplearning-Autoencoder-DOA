%%
%   Programmer: ZZF
%   Date: 2021.9.10
%
%% Multiple Rx USRPs Initialization
% close all;clc
% release(spectrumScope);
% release(radio);
radio = comm.SDRuReceiver(...
    'Platform',             'X300', ...
    'IPAddress',            '192.168.1.20,192.168.1.21,192.168.1.22,192.168.1.23', ...
    'MasterClockRate',      200e6, ...
    'CenterFrequency',      2.4e9, ...
    'Gain',                 20, ...
    'DecimationFactor',     500, ...
    'SamplesPerFrame',      4000, ...
    'OutputDataType',       'double', ...
    'ClockSource',          'External', ...
    'PPSSource',            'External', ...
    'ChannelMapping',       [1 2 3 4 5 6 7 8] );

%% Rx antenna setting
cf = 2.4e9; % 10e6 2.45e9
lambda = physconst('LightSpeed') / cf;
M0 = 8; % Total number of antenna elements
nsnapshot = radio.SamplesPerFrame; % number of snapshots

K = 2; % # of source signal
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

%% Flag setting
close all
f_normalize = 0;
f_phaseEst = 0;
closeMeasurePower = 1;

countScope = 0;
f1 = figure;
% Receiver start !
% Inside a while loop, receive the sine wave using the rx System object.
% Normalize the signal with respect to amplitude for each receive channel.
% Compute the fast fourier transform (FFT) of each normalized signal100.
% Calculate the phase difference between channel 1 and channel 2, channel 1 and channel 3, and channel 1 and channel 4.
phaseCompensatedData = zeros(nsnapshot, M0);
while 1
    [data,len] = step(radio);

    if len > 0
        if f_phaseEst
            [estimatedPhaseOffset]=phase_corr(data,f_normalize);
        end
        % phase offset compensation
        phaseCompensatedData = data;
        for i = 1:M0-1
            phaseCompensatedData(:,i+1) = estimatedPhaseOffset{i}(data(:,i+1));
        end
                timeScope(real(phaseCompensatedData));

        
        % Enable power measurement.
        if mod(countScope,100) == 0
            closeMeasurePower = 0;
            clf(f1)
        end
        %% Doing DOA
        if ~closeMeasurePower
            closeMeasurePower = 1;
            disp('DOA estimation!')
            %%ind the spacial covariance matrix,Rxx, of the received signal
            % use the sample average hat{Rxx} to estimate the Rxx
            phaseCompensatedData_r = phaseCompensatedData.';
            h_Rxx = phaseCompensatedData_r*phaseCompensatedData_r'/nsnapshot;
            save("degree_minus_40_40.mat","phaseCompensatedData_r","h_Rxx");
            DOA_run(h_Rxx,lambda,M0,K,L,L_fb,m);
        end
    end
    countScope = countScope + 1;
end
%%
% release(timeScope);
% release(spectrumScope);
% release(radio);
