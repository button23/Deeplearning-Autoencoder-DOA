%% Initialization
close all
% clear
% clc
%% Input 1
% rf_data_combine_4frames_offset = double(rf_data_combine_4frames_offset);
input_wlan = VarName4;  %% matlabdata  dspdata  VarName2 rf_data_combine_4frames_scale
input_reshape = reshape(input_wlan,2,numel(input_wlan)/2);
real_input = input_reshape(2,:).';
imag_input = input_reshape(1,:).';
RX_Data = real_input+1i*imag_input;
RX_Data_int32 = int32(RX_Data);
% rf_data_combine_4frames_arigi = input_wlan(1:1+23679);

RX_Data_int32 = dataset_wifi;
framInexx = 1;
nnn = 200;

%% Input 2
% input_wlan2 = VarName2*4/2^6;
%% Constellation Map (Data from DSP)
% frame_ca = cell(1,4);
% denom_ca = cell(1,4);
% result_ca = cell(1,4);
% result_ca_re = zeros(1924,4);
% 
% RX_Data_re = reshape(RX_Data,52,148);
% for i=1:4
% frame_ca{1,i} = RX_Data_re(:,37*(i-1)+1:37*i);
% denom_ca{1,i} = input_wlan2(52*(i-1)+1:52*i);
% result_ca{1,i} = frame_ca{1,i}./denom_ca{1,i};
% end
% for i=1:4
% result_ca_re(:,i) = reshape(result_ca{1,i},1924,1);
% end
% result_ca_re_ = int16(reshape(result_ca_re,7696,1));
% 
% figure(3)
% % xlim([-1.5 1.5])
% % ylim([-1.5 1.5])
% 
% sz = 25;
% scatter(real(result_ca_re_),imag(result_ca_re_),sz,'*')
% 
% 
% xlabel('In-phase')
% ylabel('Quadrature')
% title('Received 256 QAM Constellation')


%% DSP 32bit input data
input_wlan_int32 = abs(int32(input_wlan));
input_wlan_binary = dec2bin(input_wlan_int32);

%% wlan waveform
% cfgVHT = wlanVHTConfig('ChannelBandwidth','CBW20','MCS',8);
% % cfgVHT
% % Coded bits per single carrier for each spatial steam
% Nbpscs = 8;
% R = 3/4;
% % Coded bits per symbol (52 subcarrier * 8bits(256 QAM))
% Ncbps = 52 * Nbpscs;
% % Data bits per OFDM symbol (416 * 3/4 = 312 bits)
% Ndbps = Ncbps * R;
% % Used for scaling Matlab data in order to move to DSP
% scalingFactor = 2^13;
% % Default value is 1024
% cfgVHT.APEPLength = 1024;
% % Create random bit stream
% psdu = randi([0 1],cfgVHT.PSDULength*8,1);
% % Convert psdu into ppdu
% wholeWaveform = wlanWaveformGenerator(psdu,cfgVHT);
% wholeWaveform_3 = [wholeWaveform;wholeWaveform;wholeWaveform];

%% Create a phase and frequency offset object and introduce a 2 Hz frequency offset.
% tgacChan = wlanTGacChannel('SampleRate',20e6,'ChannelBandwidth','CBW20','DelayProfile','Model-C');
% 
% chNoise = comm.AWGNChannel('NoiseMethod','Signal to noise ratio (SNR)',...
%     'SNR',40);
% pfOffset = comm.PhaseFrequencyOffset('SampleRate',20e6,'FrequencyOffsetSource','Input port');
% % RX_Data = chNoise(tgacChan(wholeWaveform_3));
% % RX_Data = chNoise(wholeWaveform_3);
% % RX_Data = pfOffset(wholeWaveform_3,1000);
% % RX_Data = tgacChan(chNoise(pfOffset(wholeWaveform_3,2000)));
% % RX_Data = chNoise(pfOffset(wholeWaveform_3,2000));
% % RX_Data = wholeWaveform;


%% WLAN Function Synchronization
% [startOffset,meric] = wlanSymbolTimingEstimate(RX_Data,'CBW20'); %% rxSig  RX_Data WIFI_Waveform
% RX_Data_ = RX_Data(startOffset+1-length(tx)*(nnn):startOffset-length(tx)*(nnn-1));
% figure(22)
% plot(meric)

%% L-STF
% seqShort = sqrt(1/2)*[0;0;1+j;0;0;0;-1-j;0;0;0;1+j;0;0;0;-1-j;
%     0;0;0;-1-j;0;0;0;1+j;0;0;0;0;0;0;0;-1-j;
%     0;0;0;-1-j;0;0;0;1+j;0;0;0;1+j;0;0;0;1+j;
%     0;0;0;1+j;0;0];
% STFpriorIFFT = [zeros(6,1);seqShort;zeros(5,1)];
% STF_ifft = ifft(ifftshift(STFpriorIFFT,1))*(1/sqrt(12))*64;
% STFifft_double = [STF_ifft;STF_ifft];
% STFifft_full = [STFifft_double(97:128);STFifft_double];
%
% % [r,lags] = xcorr(input,STF_ifft(1:16));
% [r,lags] = xcorr(RX_Data,STFifft_double);
%
% % figure(12)
% % plot(abs(r((length(RX_Data)-length(STF_ifft)-1):end))/max(abs(r)))
%
% % figure(10)
% % plot(real(STF_ifft));
% % xlabel('Samples')
% % ylabel('Real')
% % figure(11)
% % plot(imag(STF_ifft));
% % xlabel('Samples')
% % ylabel('Imaginary')
%% L-LTF
seqLong = [1;1;-1;-1;1;1;-1;1;-1;1;1;1;1;1;1;
    -1;-1;1;1;-1;1;-1;1;1;1;1;0;1;-1;-1;1;
    1;-1;1;-1;1;-1;-1;-1;-1;-1;1;1;-1;-1;1;-1;
    1;-1;1;1;1;1];
LTFpriorIFFT = [zeros(6,1);seqLong;zeros(5,1)];
LTF_ifft = ifft(ifftshift(LTFpriorIFFT,1))*(1/sqrt(52))*64;
LTFifft_double = [LTF_ifft;LTF_ifft];
LTFifft_full = [LTFifft_double(97:128);LTFifft_double];

% [r,lags] = xcorr(RX_Data,LTFifft_full);
% figure(14)
% plot(abs(RX_Data(2960*30:2960*50)))  %%rxSig RX_Data

% data_test = RX_Data(817:817+2959);
% data_test(641:720) = data_test(641:720) * (64/sqrt(52));

% % r_corr = r((length(RX_Data)-length(LTFifft_full)+1):end)/max(abs(r)); % pick the values aside from zeros
% % figure(13)
% % plot(abs(r_corr))
% 
% % index = find(abs(r_corr)>0.8);
% % for ii=1:length(index)
% %     if abs(r_corr(index(ii)+64)) > 0.8  % make sure it is the first peak, not the truncted frame
% %         offset_ltf = index(ii) - 256;  % the beginning of LTF
% %         break
% %     end
% % end
% 
% 
% r_corr = r(5920:5920+2960+160-1)/max(abs(r)); % pick the values aside from zeros
% figure(13)
% plot(abs(r_corr))
% index = find(abs(r_corr) == max(abs(r_corr)));
% offset_ltf = index - 160;  % STF 160 sample, LTF 160 sample
% 
% frame_sync = RX_Data(offset_ltf:offset_ltf+2959);
% 
% % frame_sync_offset = pfOffset(frame_sync,4000);

%% DSP coding (for synchronization)
scaling = 2^5;
synchResult = zeros(2960+160+1,1);
% synchResult = zeros(5040+160+1,1);

% % header LTF sequence (after ifft)
% realLTF_int16 = int16(real(LTFifft_full) * scaling);
% imagLTF_int16 = int16(imag(LTFifft_full) * scaling);
% head_ltf = [imagLTF_int16.';realLTF_int16.'];
% head_ltf_reshape = reshape(head_ltf,320,1);


realLTF = real(LTFifft_full) * scaling;
imagLTF = imag(LTFifft_full) * scaling;

for dd=1:length(synchResult) %% only need to calculate 2960 times which is half of the two frames.
    mid = RX_Data(dd:dd+160-1);  %%synchFrame_fourframe_offset
    rrrMid = real(mid);
    iiiMid = imag(mid);
    realll = 0;
    imaggg = 0;
    for ccc = 1:160
        realll =  realll +(rrrMid(ccc) * realLTF(ccc) + iiiMid(ccc) * imagLTF(ccc));
        imaggg =  imaggg +(rrrMid(ccc) * imagLTF(ccc) - realLTF(ccc) * iiiMid(ccc));
    end
    
    synchResult(dd) = abs((realll)) + abs((imaggg));
end
% figure(16)
% plot((synchResult))

% index_ = find(abs(synchResult) == max(abs(synchResult)));

maxx = 0;
index_c = 1;
for nn = 1:length(synchResult)
    if abs(synchResult(nn)) > maxx
        maxx = abs(synchResult(nn));
        index_c = nn;
    end
end

% synchResult_offset = index_c - 160
% if synchResult_offset < 0 
%     synchResult_offset = index_c + 4880 
      synchResult_offset = index_c + 2800 
% end

% synchFrame = RX_Data((synchResult_offset + (framInexx-1) * 5040):(synchResult_offset+5039 + (framInexx-1) * 5040)); %% the same as  frame_sync which has different noise in it
synchFrame = RX_Data((synchResult_offset + (framInexx-1) * 2960):(synchResult_offset+2959 + (framInexx-1) * 2960)); %% the same as  frame_sync which has different noise in it


figure(17)
plot(abs(synchFrame))


%% Waveform Figure
% x = 0:0.05:148-0.005;
% figure(52)
% plot(x,(frame_sync))
% xlabel('Time (us)')
% ylabel('Amplitude')
% title('Received Signal (Time Domain)')
% xlim([0 148])
%
%
% figure(53)
% [Pxx_wifi,F_wifi] = periodogram(frame_sync,[],length(frame_sync),20000000,'centered'); %wholeWaveform  frame_sync
% plot(F_wifi/10^6,10*log10(Pxx_wifi))
% xlabel('Frequency (MHz)')
% ylabel('Magnitude (dB)')
% title('Received Signal (Frequency Domain)')
% xlim([-12 12])


%% PSD Figure
% Fs = 200000; % sampling rate
% N = length(frame_sync);
% frame_sync_reshape = reshape(frame_sync,80,37);
% frame_no_cp = frame_sync_reshape(17:end,:);
% frame_single = frame_no_cp(:,35);
% frame_single_autocorr = xcorr(frame_single);
% frame_single_autocorr_fft = fftshift(fft(frame_single_autocorr));
% figure(55)
% plot(10*log10(abs(frame_single_autocorr_fft)))
% xlabel('Frequency (MHz)')
% ylabel('Magnitude (dB)')
% title('Received Signal (Frequency Domain)')
%  
%% PSD Figure
% x = 0:63; % FFT point
% xx = x*0.3125;  % The frequency of first harmonic wave is 312.5kHz %0.3125 means the unit of x coordinate is MHz
% wifi_autocorr = xcorr(frame_sync);
% wifi_autocorr_fft = fftshift(fft(wifi_autocorr));
% figure(54)
% plot(xx,10*log10(abs(wifi_autocorr_fft)))
% xlabel('Frequency (MHz)')
% ylabel('Magnitude (dB)')
% title('Received Signal (Frequency Domain)')



%% Preamble
% preamble = [STFifft_full;LTFifft_full];
%
% PREAMBLEcorr = xcorr(preamble,LTF_ifft);
% figure(9)
% plot(abs(PREAMBLEcorr(257:end))/max(abs(PREAMBLEcorr)));

%% Timing Synchronization
% [startOffset,metrrr] = wlanSymbolTimingEstimate(RX_Data,'CBW20');
% sync_frame = RX_Data(startOffset+1:startOffset+2960);
% % figure(23)
% % plot(metrrr)
%
% isequal(startOffset,offset_ltf)
