%%
%   Programmer: ZZF
%   Date: 2021.9.10
%
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
    'ChannelMapping',       [1]);

%% Transmit start !
% Generate a sine wave of 30 kHz for transmission.
% The sample rate is calculated from the master clock rate
% and interpolation factor specified for an X310 radio System object configuration.
% Set the output data type of the sine wave as 'double'.

sinewave = dsp.SineWave(1,1e3);
sinewave.SampleRate = radio.MasterClockRate/radio.InterpolationFactor; % Ns = Fs/Fw=(200e6/500)/1e3=400
sinewave.SamplesPerFrame = 4000;
sinewave.OutputDataType = 'double';
sinewave.ComplexOutput = true;
data = step(sinewave);

% Set the frame duration for the sine wave to transmit based on the samples per frame and sample rate.
% Create time scope and frequency scope System objects to display time-domain and frequency-domain signals,
% respectively. Display a message when transmission starts.

frameDuration = (sinewave.SamplesPerFrame)/(sinewave.SampleRate);
time = 0;
% timeScope = timescope('SampleRate',sinewave.SampleRate);
timeScope = timescope('TimeSpanSource','Property','TimeSpan',frameDuration,...
    'SampleRate',sinewave.SampleRate);
% spectrumScope = dsp.SpectrumAnalyzer('SampleRate',sinewave.SampleRate);
disp("Transmission Started");
timeScope((data));

% spectrumScope(ltee.waveform);
% %% weight multiplication
% % clear all;clc;
% c = physconst('LightSpeed');        % propagation speed
% fc = 2.4e9;      % carrier frequency
% lambda = c/fc;  % wavelength
% antenna = phased.IsotropicAntennaElement(...
%     'FrequencyRange',[3e8 3e9]);
% txarray = phased.ULA('NumEldataements',2,'ElementSpacing',lambda/2);
% % viewArray(txarray);
% txmipos_realPos = getElementPosition(txarray);
% txmipos = txmipos_realPos/lambda;  % x, y, z, the default coordinate that used is y
% % the reason divide by lambda is normalized
% %%
% %
% txarraystv = phased.SteeringVector('SensorArray',txarray,'PropagationSpeed',c);
% txang = [30;0];  % tx angle (azimuth; elevation)          -180<=azimuth<=180  -90<=elevation<=90
% wt = txarraystv(fc,txang)';
% wt = wt.'
% 
% %% lte data
% data = ltee.waveform;
% % data_long = data;
% % data_long = []
% % for aaa = 1 : 10
% %     data_long = [data_long;data];
% % end
% 
% %% test
% %% Cell-Wide Settings
% %% Cell-Wide Settings
% rmccfg.RC = 'R.7';
% ncodewords = 1;
% enb = lteRMCDL(rmccfg,ncodewords);
% enb.NSubframe = 0;
% cellRsInd = lteCellRSIndices(enb);
% 
% [offset,corr] = lteDLFrameOffset(enb,data); %%txWaveform rxWaveform waveform
% figure
% plot(corr)
% %% weight multiplication
% data_raw_1 = data.*wt(1);
% data_raw_2 = data.*wt(2);
% send_data = [data_raw_1, data_raw_2];
% send_data_scale = send_data*2^14;
% disp([max(real(send_data_scale)) max(imag(send_data_scale))])
%
% Inside a while loop, transmit the sine wave using the tx System object.
% Display a message when transmission is complete. Release the radio System object.
while true
    radio(data);  %data  send_data_scale
    %     time = time+frameDuration;
end
%
% disp("Transmission Stopped");
