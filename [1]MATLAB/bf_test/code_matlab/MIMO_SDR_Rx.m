%%
%   Programmer: ZZF
%   Date: 2021.9.10
% 
%% Connect USRPs
% connectedRadios = findsdru;
% %%
% n_usrp = length(connectedRadios);
% address = cell(n_usrp,1);
% platform = cell(n_usrp,1);
% 
% for i = 1 : n_usrp
%     if strncmp(connectedRadios(i).Status, 'Success', 7)
% %         address{i} = connectedRadios(i).IPAddress;
%         address{i} = '192.168.1.23';
% %         platform{i} = connectedRadios(i).Platform;
%         platform{i} = 'X300';
% 
%     else
%         disp("Device is not found !!")
%         quit
%     end
% end

%% Initialization
prmQPSKReceiver = receiver_init('X300');
prmQPSKReceiver.Platform = 'X300';
prmQPSKReceiver.Address = '192.168.1.23';

%% Receiver start !
% qpskRx = QPSKReceiver(...
%     'ModulationOrder',                      prmQPSKReceiver.ModulationOrder, ...
%     'SampleRate',                           prmQPSKReceiver.Fs, ...
%     'DecimationFactor',                     prmQPSKReceiver.Decimation, ...
%     'FrameSize',                            prmQPSKReceiver.FrameSize, ...
%     'HeaderLength',                         prmQPSKReceiver.HeaderLength, ...
%     'NumberOfMessage',                      prmQPSKReceiver.NumberOfMessage, ...
%     'PayloadLength',                        prmQPSKReceiver.PayloadLength, ...
%     'DesiredPower',                         prmQPSKReceiver.DesiredPower, ...
%     'AveragingLength',                      prmQPSKReceiver.AveragingLength, ...
%     'MaxPowerGain',                         prmQPSKReceiver.MaxPowerGain, ...
%     'RolloffFactor',                        prmQPSKReceiver.RolloffFactor, ...
%     'RaisedCosineFilterSpan',               prmQPSKReceiver.RaisedCosineFilterSpan, ...
%     'InputSamplesPerSymbol',                prmQPSKReceiver.Interpolation, ...
%     'MaximumFrequencyOffset',               prmQPSKReceiver.MaximumFrequencyOffset, ...
%     'PostFilterOversampling',               prmQPSKReceiver.Interpolation/prmQPSKReceiver.Decimation, ...
%     'PhaseRecoveryLoopBandwidth',           prmQPSKReceiver.PhaseRecoveryLoopBandwidth, ...
%     'PhaseRecoveryDampingFactor',           prmQPSKReceiver.PhaseRecoveryDampingFactor, ...
%     'TimingRecoveryDampingFactor',          prmQPSKReceiver.TimingRecoveryDampingFactor, ...
%     'TimingRecoveryLoopBandwidth',          prmQPSKReceiver.TimingRecoveryLoopBandwidth, ...
%     'TimingErrorDetectorGain',              prmQPSKReceiver.TimingErrorDetectorGain, ...
%     'PreambleDetectorThreshold',            prmQPSKReceiver.PreambleDetectorThreshold, ...    
%     'DescramblerBase',                      prmQPSKReceiver.ScramblerBase, ...
%     'DescramblerPolynomial',                prmQPSKReceiver.ScramblerPolynomial, ...
%     'DescramblerInitialConditions',         prmQPSKReceiver.ScramblerInitialConditions,...
%     'BerMask',                              prmQPSKReceiver.BerMask);

%    
radio = comm.SDRuReceiver(...
    'Platform',             prmQPSKReceiver.Platform, ...
    'IPAddress',            prmQPSKReceiver.Address, ...
    'MasterClockRate',      prmQPSKReceiver.MasterClockRate, ...
    'CenterFrequency',      prmQPSKReceiver.USRPCenterFrequency, ...
    'Gain',                 prmQPSKReceiver.USRPGain, ...
    'DecimationFactor',     prmQPSKReceiver.USRPDecimationFactor, ...
    'SamplesPerFrame',      prmQPSKReceiver.USRPFrameLength, ...
    'OutputDataType',       'double', ...
    'ClockSource',          prmQPSKReceiver.Clock, ...
    'PPSSource',            prmQPSKReceiver.PPS, ...
    'ChannelMapping',       [1 2] );
      
% %%
% % Initialize variables
% len = uint32(0);
% rcvdSignal = complex(zeros(prmQPSKReceiver.USRPFrameLength,1));
% BER = [];
% currentTime = 0;
% 
% while currentTime < prmQPSKReceiver.StopTime
%     % Keep accessing the SDRu System object output until it is valid
%     while len <= 0
%         [rcvdSignal, len] = step(radio);
%     end
%     
%     % When the SDRu System object output is valid, decode the received
%     % message
%     if len > 0
%         [~, ~, ~, BER] = qpskRx(rcvdSignal); % Receiver
%     end
%     
%     len=uint32(0);
%     currentTime = currentTime + prmQPSKReceiver.USRPFrameTime;
% end
% 
% scatterplot(rcvdSignal)
% 
% 
% % timeScope = timescope('TimeSpanSource','Property','TimeSpan',4/30e3,...
% %                       'SampleRate',100e6/100);
% % spectrumScope = dsp.SpectrumAnalyzer('SampleRate',prmQPSKReceiver.USRPFrontEndSampleRate); 
% % % disp("Transmission Started"); 
% % timeScope(rcvdSignal); 
% % spectrumScope(rcvdSignal);

%% Phase check
frameduration = 1e4/(200e6/500); 
time = 0; 
timeScope = timescope('TimeSpanSource','Property','TimeSpan',...
                      5/30e3,'SampleRate',200e6/200);
spectrumScope = dsp.SpectrumAnalyzer('SampleRate',200e6/500); 
spectrumScope.ReducePlotRate = true; 
disp("Reception Started");

% Inside a while loop, receive the sine wave using the rx System object. 
% Normalize the signal with respect to amplitude for each receive channel. 
% Compute the fast fourier transform (FFT) of each normalized signal. 
% Calculate the phase difference between channel 1 and channel 2, channel 1 and channel 3, and channel 1 and channel 4.

while time < 10000 
 [data,len] = step(radio); 
 if len > 0 
    amp(1) = max(abs(data(:,1))); 
    amp(2) = max(abs(data(:,2))); 
%     amp(3) = max(abs(data(:,3))); 
%     amp(4) = max(abs(data(:,4))); 
    maxAmp = max(amp); 
    if any(~amp)  
       NormalizedData = data; 
    else 
      NormalizedData(:,1) = maxAmp/amp(1)*data(:,1); 
      NormalizedData(:,2) = maxAmp/amp(2)*data(:,2); 
%       NormalizedData(:,3) = maxAmp/amp(3)*data(:,3); 
%       NormalizedData(:,4) = maxAmp/amp(4)*data(:,4); 
    end 
    freqOfFirst = fft(NormalizedData(:,1)); 
    freqOfSecond = fft(NormalizedData(:,2)); 
%     freqOfThird = fft(NormalizedData(:,3)); 
%     freqOfFourth = fft(NormalizedData(:,4)); 
    angle1 = rad2deg(angle(max(freqOfFirst)/max(freqOfSecond))); 

%     % phase offset compensation
%     pfo = comm.PhaseFrequencyOffset('PhaseOffset',angle1);
%     NormalizedData(:,2) = pfo(NormalizedData(:,2));

        % test 
%     kkk = NormalizedData(:,1);
%     ttt = NormalizedData(:,2);
%     freqofsss = fft(ttt);
%     angle2 = rad2deg(angle(max(freqOfFirst)/max(freqofsss))); 

    
%     angle2 = rad2deg(angle(max(freqOfFirst)/max(freqOfSecond))); 

%     if time%250 == 0
%         disp(angle1)
%     end
%     angle2 = rad2deg(angle(max(freqOfFirst)/max(freqOfThird))); 
%     angle3 = rad2deg(angle(max(freqOfFirst)/max(freqOfFourth))); 
%     timeScope([real(NormalizedData),imag(NormalizedData)]); 
    timeScope([real(NormalizedData(:,1)),real(NormalizedData(:,2))]); 

    spectrumScope(NormalizedData); 
 end 
    time = time + frameduration; 
end

% Display the calculated phase difference between channel 1 and each of the other channels of TwinRX daughterboard.

disp([' Phase difference between channel 1 and 2: ', num2str(angle1)]); 
% disp([' Phase difference between channel 1 and 3: ', num2str(angle2)]); 
% disp([' Phase difference between channel 1 and 4: ', num2str(angle3)]); 
disp("Reception Ended"); 
release(timeScope); 
release(spectrumScope); 
release(radio);

% release(qpskRx);
% release(radio);
