% Generate a sine wave of 30 kHz for transmission. 
% The sample rate is calculated from the master clock rate 
% and interpolation factor specified for an N210 radio System object configuration. 
% Set the output data type of the sine wave as 'double'.

sinewave = dsp.SineWave(1,30e3); 
sinewave.SampleRate = 200e6/50; 
sinewave.SamplesPerFrame = 50e4; 
sinewave.OutputDataType = 'double'; 
sinewave.ComplexOutput = true;

data = step(sinewave);

% Set the frame duration for the sine wave to transmit based on the samples per frame and sample rate. 
% Create time scope and frequency scope System objects to display time-domain and frequency-domain signals, 
% respectively. Display a message when transmission starts.

frameDuration = (sinewave.SamplesPerFrame)/(sinewave.SampleRate); 
time = 0;
timeScope = timescope('TimeSpanSource','Property','TimeSpan',10/30e3,...
                      'SampleRate',100e6/100);
spectrumScope = dsp.SpectrumAnalyzer('SampleRate',sinewave.SampleRate); 
disp("Transmission Started"); 
timeScope(data); 
spectrumScope(data);
%%
% Inside a while loop, transmit the sine wave using the tx System object. Display a message when transmission is complete. Release the radio System object.

while time<30
    tx(data); 
    time = time+frameDuration; 
end 
disp("Transmission Stopped"); 
release(tx);