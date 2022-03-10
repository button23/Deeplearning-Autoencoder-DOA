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

sinewave = dsp.SineWave(10,1e3);
sinewave.SampleRate = 200e6/500; % Ns = Fs/Fw=(200e6/500)/1e3=400
sinewave.SamplesPerFrame = 4000;
sinewave.OutputDataType = 'double';
sinewave.ComplexOutput = true;
data = step(sinewave);
% data_two = [data,data];

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

% Inside a while loop, transmit the sine wave using the tx System object.
% Display a message when transmission is complete. Release the radio System object.
while true
    radio(data);  %data  send_data_scale data_two
    %     time = time+frameDuration;
end
%
% disp("Transmission Stopped");
