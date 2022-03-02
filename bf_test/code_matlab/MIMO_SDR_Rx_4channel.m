%%
%   Programmer: ZZF
%   Date: 2021.9.10
%
%% Multiple Rx USRPs Initialization
close all;clc
% release(spectrumScope);
% release(radio);
radio = comm.SDRuReceiver(...
    'Platform',             'X300', ...
    'IPAddress',            '192.168.1.20,192.168.1.21,192.168.1.22,192.168.1.23', ...
    'MasterClockRate',      200e6, ...
    'CenterFrequency',      2.4e9, ...
    'Gain',                 23, ...
    'DecimationFactor',     100, ...
    'SamplesPerFrame',      4000, ...
    'OutputDataType',       'double', ...
    'ClockSource',          'External', ...
    'PPSSource',            'External', ...
    'ChannelMapping',       [1 2 3 4 5 6 7 8] );


%% Tx antenna setting
cf = 2.4e9; % 10e6 2.45e9
lambda = physconst('LightSpeed') / cf;
ula = phased.ULA('NumElements',2,'ElementSpacing',lambda/2);

% %% MVDR beamformer (also called the Capon beamformer)
% % The MVDR beamformer preserves the signal arriving along a desired
% % direction, while trying to suppress signals coming from other directions.
% % 'TrainingInputPort' means we have access to the interence data
% mvdrbeamformer = phased.MVDRBeamformer('SensorArray',ula,...
%     'OperatingFrequency',cf,'Direction',angle,...
%     'WeightsOutputPort',true,'TrainingInputPort',true);
% [yMVDR,wMVDR] = mvdrbeamformer(rx_signal,rx_int);
%
% % Looking at the response pattern of the beamformer, we see two deep nulls
% % along the interference directions, (30 and 50 degrees).
% % The beamformer also has a gain of 0 dB along the target direction
% % of 45 degrees. Thus, the MVDR beamformer preserves the target signal
% % and suppresses the interference signals.  PhaseShift pattern does not
% % null out the interference at al
% figure
% plot(t,abs(yMVDR)); axis tight;
% title('Output of MVDR Beamformer With Presence of Interference');
% xlabel('Time (s)');ylabel('Magnitude (V)');
%
% figure
% % Plot array response with phase-shft weighting for each antenna element
% pattern(ula,cf,-180:180,0,'Weights',wMVDR,'Type','powerdb',...
%     'PropagationSpeed',physconst('LightSpeed'),'Normalize',false,...
%     'CoordinateSystem','rectangular')
% axis([-90 90 -80 20])
%
% hold on;
% pattern(ula,cf,-180:180,0,'Weights',wPS,'Type','powerdb',...
%     'PropagationSpeed',physconst('LightSpeed'),'Normalize',false,...
%     'CoordinateSystem','rectangular')
% hold off;
% legend('MVDR','Phase Shift')

%% Receiver start !
% release(timeScope);

frameduration = radio.SamplesPerFrame/(radio.MasterClockRate/radio.DecimationFactor);
time = 0;
timeScope = timescope('TimeSpanSource','Property','TimeSpan',...
    4000/1e3,'SampleRate',radio.MasterClockRate/radio.DecimationFactor);
% spectrumScope = dsp.SpectrumAnalyzer('SampleRate',200e6/500);
% spectrumScope.ReducePlotRate = true;
disp("Reception Started");
closePhaseEstimate = 1;
closeMeasurePower = 1;
closeNormalize = 1;

flagPowerMeasure = 0;
flagPhase = 0;
countScope = 0;
% figure

f3 = figure;

%%Inside a while loop, receive the sine wave using the rx System object.
% Normalize the signal with respect to amplitude for each receive channel.
% Compute the fast fourier transform (FFT) of each normalized signal.
% Calculate the phase difference between channel 1 and channel 2, channel 1 and channel 3, and channel 1 and channel 4.
while time < 1000000
    [data,len] = step(radio);
    if len > 0
        amp(1) = max(abs(data(:,1)));
        amp(2) = max(abs(data(:,2)));
        amp(3) = max(abs(data(:,3)));
        amp(4) = max(abs(data(:,4)));
        amp(5) = max(abs(data(:,5)));
        amp(6) = max(abs(data(:,6)));
        amp(7) = max(abs(data(:,7)));
        amp(8) = max(abs(data(:,8)));

        NormalizedData = data;
        % Normalize the data to make the data of all the channel have the
        % same amplitude.
        if ~closeNormalize
            maxAmp = max(amp);
            NormalizedData(:,1) = maxAmp/amp(1)*data(:,1);
            NormalizedData(:,2) = maxAmp/amp(2)*data(:,2);
            NormalizedData(:,3) = maxAmp/amp(3)*data(:,3);
            NormalizedData(:,4) = maxAmp/amp(4)*data(:,4);
            NormalizedData(:,5) = maxAmp/amp(5)*data(:,5);
            NormalizedData(:,6) = maxAmp/amp(6)*data(:,6);
            NormalizedData(:,7) = maxAmp/amp(7)*data(:,7);
            NormalizedData(:,8) = maxAmp/amp(8)*data(:,8);
        end
        % Estimate the phase difference between channels.
        if ~closePhaseEstimate
            if ~flagPhase
                freq1 = fft(NormalizedData(:,1));
                freq2 = fft(NormalizedData(:,2));
                freq3 = fft(NormalizedData(:,3));
                freq4 = fft(NormalizedData(:,4));
                freq5 = fft(NormalizedData(:,5));
                freq6 = fft(NormalizedData(:,6));
                freq7 = fft(NormalizedData(:,7));
                freq8 = fft(NormalizedData(:,8));

                angle1 = rad2deg(angle(max(freq1)/max(freq2)));
                angle2 = rad2deg(angle(max(freq1)/max(freq3)));
                angle3 = rad2deg(angle(max(freq1)/max(freq4)));
                angle4 = rad2deg(angle(max(freq1)/max(freq5)));
                angle5 = rad2deg(angle(max(freq1)/max(freq6)));
                angle6 = rad2deg(angle(max(freq1)/max(freq7)));
                angle7 = rad2deg(angle(max(freq1)/max(freq8)));

                % phase offset compensation
                pfo1 = comm.PhaseFrequencyOffset('PhaseOffset',angle1);
                pfo2 = comm.PhaseFrequencyOffset('PhaseOffset',angle2);
                pfo3 = comm.PhaseFrequencyOffset('PhaseOffset',angle3);
                pfo4 = comm.PhaseFrequencyOffset('PhaseOffset',angle4);
                pfo5 = comm.PhaseFrequencyOffset('PhaseOffset',angle5);
                pfo6 = comm.PhaseFrequencyOffset('PhaseOffset',angle6);
                pfo7 = comm.PhaseFrequencyOffset('PhaseOffset',angle7);
                flagPhase = 1;
                disp('Phase corrected !')
            end
        end
        NormalizedData(:,2) = pfo1(NormalizedData(:,2));
        NormalizedData(:,3) = pfo2(NormalizedData(:,3));
        NormalizedData(:,4) = pfo3(NormalizedData(:,4));
        NormalizedData(:,5) = pfo4(NormalizedData(:,5));
        NormalizedData(:,6) = pfo5(NormalizedData(:,6));
        NormalizedData(:,7) = pfo6(NormalizedData(:,7));
        NormalizedData(:,8) = pfo7(NormalizedData(:,8));

        % Update the timescope.
        if mod(countScope,500) == 0
            timeScope(real(NormalizedData));
            disp('scope update !')
%             flagPhase = 0;
        end
        % Enable phase estimation.
        if ~closePhaseEstimate
            if mod(countScope,2000) == 0
                flagPhase = 0;
            end
        end

        if mod(countScope,100) == 0
            closeMeasurePower = 0;
            clf(f3)
        end
%         Enable power measurement.

        if ~closeMeasurePower
            closeMeasurePower = 1;
            disp('DOA estimation!')
            
                                        %% Find the spacial covariance matrix,Rxx, of the received signal
                            % use the sample average hat{Rxx} to estimate the Rxx
                            NormalizedData_ = NormalizedData.';
                            h_Rxx = NormalizedData_*NormalizedData_'/4000;
            
        end

        %     spectrumScope(NormalizedData);
    end
    time = time + frameduration;
    countScope = countScope + 1;
end
%

%%
% release(timeScope);
% release(spectrumScope);
% release(radio);
