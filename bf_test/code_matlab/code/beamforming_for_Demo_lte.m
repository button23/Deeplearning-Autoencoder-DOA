%%
%   Programmer: ZZF
%   Date: 2021.9.10
%
%% Multiple Rx USRPs Initialization
close all;clc
% release(timeScope);
% release(spectrumScope);
% release(radio);
radio = comm.SDRuReceiver(...
    'Platform',             'X300', ...
    'IPAddress',            '192.168.1.24', ...
    'MasterClockRate',      184.32e6, ...
    'CenterFrequency',      2.4e9, ...
    'Gain',                 10, ...
    'DecimationFactor',     12, ...
    'SamplesPerFrame',      15360*10, ...
    'ClockSource',          'External', ...
    'PPSSource',            'External', ...
    'ChannelMapping',       [1] );
%% Channel Estimator Configuration
cec.PilotAverage = 'UserDefined'; % Pilot averaging method
cec.FreqWindow = 9;               % Frequency averaging window in REs
cec.TimeWindow = 9;               % Time averaging window in REs
cec.InterpType = 'Cubic';         % Cubic interpolation
cec.InterpWinSize = 3;            % Interpolate up to 3 subframes
% simultaneously
cec.InterpWindow = 'Centred';     % Interpolation windowing method
%% Cell-Wide Settings
rmccfg.RC = 'R.7';
ncodewords = 1;
enb = lteRMCDL(rmccfg,ncodewords);
enb.NSubframe = 0;
cellRsInd = lteCellRSIndices(enb);

%% Receiver start !
demoFunc = 0;
frameduration = 1;
time = 0;
timeScope = timescope('TimeSpanSource','Property','TimeSpan',...
    frameduration,'SampleRate',184.32e6/12);
spectrumScope = dsp.SpectrumAnalyzer('SampleRate',184.32e6/12);
% spectrumScope.ReducePlotRate = true;
disp("Reception Started");
num_iteration = 20;
% figure
% Inside a while loop, receive the sine wave using the rx System object.
% Normalize the signal with respect to amplitude for each receive channel.
% Compute the fast fourier transform (FFT) of each normalized signal.
% Calculate the phase difference between channel 1 and channel 2, channel 1 and channel 3, and channel 1 and channel 4.

while time < 100000
    [data,len] = step(radio);
    %     d_scale = 2^-4;
    data = double(data);
    
    data = imag(data) + 1i*real(data);
    
    
    %% CFO Estimation and Compensation
    cfo_offset = lteFrequencyOffset(enb,data); %waveform waveform
%         if mod(time,30) == 0
%               disp(['Frequency offset: ' num2str(cfo_offset) ' Hz'])
%         end
    waveform1_nocfo = lteFrequencyCorrect(enb,data,cfo_offset);
%         cfo_o = lteFrequencyOffset(enb,waveform1_nocfo); %waveform waveform
%             disp(['Frequency offset: ' num2str(cfo_o) ' Hz'])
    
    %% Synchronization
    [offset,corr] = lteDLFrameOffset(enb,waveform1_nocfo); %%txWaveform rxWaveform waveform
%     figure
%     plot(corr)
    fft_size = 1024; cp = 72; one_sym = fft_size + cp;
    if (offset + 15360 <= 15360*10 && offset>0)
        subframe_data = waveform1_nocfo(offset+1 :offset + 15360);
        spectrumScope(subframe_data);
        
        timeScope(abs(subframe_data));
        %         if demoFunc == 1
        %             %% OFDM Demodulation (function)
        %             rxGrid = lteOFDMDemodulate(enb,subframe_data);
        %         else
        %% OFDM Demodulation (coding)
        rxGrid_fft=reshape(subframe_data(1:15344),one_sym,14);
        rxGrid_fft_nocp = fft(rxGrid_fft(73:end,:));
        rxGrid = [rxGrid_fft_nocp(end-299:end,:);rxGrid_fft_nocp(2:301,:)];
        %         end
        %% Get Cell refrence signal
        crs_subfram = rxGrid(cellRsInd); % get crs from received signal
        
        %% Channel Estimation
        [estChannel, noiseEst] = lteDLChannelEstimate(enb,cec,rxGrid); %%fftc_output_grid  rxGrid
        
        %% MMSE Equalization
        eqGrid = lteEqualizeZF(rxGrid, estChannel);
        
        %% plot
        if mod(time,num_iteration) == 0
            figure(10)
            subplot(2,2,1)
            rxGrid_one = reshape(rxGrid,8400,1);
            rxGrid_one_nor = rxGrid_one./max(abs(rxGrid_one));
            scatter(real(rxGrid_one_nor),imag(rxGrid_one_nor));
            title("Constelation Map before Channel Compensation (64QAM)")
            xlabel("In-phase")
            ylabel("Quadrature")
            
            subplot(2,2,2)
            eqGrid_one = reshape(eqGrid,600*14,1);
            scatter(real(eqGrid_one),imag(eqGrid_one))
            xlim([-1.5 1.5])
            ylim([-1.5 1.5])
            title("Constelation Map after Channel Compensation (64QAM)")
            xlabel("In-phase")
            ylabel("Quadrature")
            
            subplot(2,2,3)
            scatter(real(crs_subfram),imag(crs_subfram))
            %             xlim([-1.5 1.5])
            %             ylim([-1.5 1.5])
            title("Reference Signal")
            xlabel("In-phase")
            ylabel("Quadrature")
            
            subplot(2,2,4)
            surf(abs(estChannel(:,:,1,1)))  %%channel_status_grid estChannel
            title("Channel Estimation")
            xlabel("Time")
            ylabel("Frequency")
            
        end
        time = time + 1;
    end
end
%%
release(timeScope);
% release(spectrumScope);
release(radio);

%% power
% powerx = 0:10:90;
% powery = [65.5,76.36,62.5,71.8,65.9,69.5,74.75,74.3,65.6,65.4];
% %reference
% power_refe = [2,1.9261,1.7182,1.4142,1.0642,0.7186,0.4178,0.1892,0.0477,0];
% figure
% plot(powerx,powery)
% hold on
% plot(powerx,power_refe)
% %% 0 degree
% powerx = 0:10:90;
% powerr = [85.5 84.8 85.5 84.8 81.6 82.7 81 77 80.7 74.8];
% powerl = [85.1 83.8 84.6 85.1 82.4 83.2 84.6 83.6 83 81.1];
% %reference
% power_refe = [2,1.9261,1.7182,1.4142,1.0642,0.7186,0.4178,0.1892,0.0477,0];
% figure
% % subplot(2,1,1)
% plot(powerx,powerr)
% % subplot(2,1,2)
% hold on
% plot(powerx,powerl)
% % hold on
% 
% %% 30 degree
% powerx = [-90,-60, -45,-30,0,0,30,45,60,90];
% powerr_30 = [83.16, 84.23, 78.417, 76.677,71.623];
% powerl_30 = [84.837,74.137,79.236,79.699,81.347];
% powerl_30_fl = fliplr(powerl_30);
% power_30 = [powerl_30_fl,powerr_30];
% figure
% plot(powerx,power_30)
% % %reference
% % power_refe = [2,1.9261,1.7182,1.4142,1.0642,0.7186,0.4178,0.1892,0.0477,0];
% % figure
% % % subplot(2,1,1)
% % plot(powerx,powerr_30)
% % % subplot(2,1,2)
% % hold on
% % plot(powerx,powerl_30)
% % % hold on
% % figure
% % plot(powerr_30)
% % hold on
% % plot(powerl_30)
% %% 60 degree
% powerx = [-90,-60, -45,-30,0,0,30,45,60,90];
% powerr_60 = [78 83.5 79.2 80.2 79.3];
% powerl_60 = [80 75 80 82.5 80.6];
% powerr_60_lr = flip(powerl_60);
% 
% power_60 = [powerr_60_lr,powerr_60];
% % power_60 = [powerl_60,powerr_60];
% 
% figure
% plot(powerx,power_60)
% %% 45 degree
% % powerx = [-90,-60, -45,-30,0,0,30,45,60,90];
% powerr_45 = [84.344,76.357,74.107,77.139,77.031];
% powerl_45 = [84.657,75.681,73.667,71.235,70.341];
% powerr_45_lr = fliplr(powerl_45);
% 
% power_45 = [powerr_45_lr,powerr_45];
% figure
% plot(powerx,power_45)

