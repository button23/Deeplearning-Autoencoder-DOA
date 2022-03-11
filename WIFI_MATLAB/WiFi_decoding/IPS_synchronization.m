%% Initialization
close all
% clear
clc

%% RF data input format conversion
input_wlan = VarName2;  %% matlabdata  dspdata  VarName2 rf_data_combine_4frames_scale
input_reshape = reshape(input_wlan,2,numel(input_wlan)/2);

% Imaginary and Real part separation
imag_input = input_reshape(1,:).';
real_input = input_reshape(2,:).';

RX_Data = real_input+1i*imag_input;
RX_Data_int32 = int32(dataset_wifi);

% choose picked frame index
nnn = 5;
sampleNum =4000;  %2960

%% L-LTF generation
% Also confirmed from matlab provided code that need to do fft and cp
% adding to the sequence. Eventually, we will have 160samples.
seqLong = [1;1;-1;-1;1;1;-1;1;-1;1;1;1;1;1;1;
    -1;-1;1;1;-1;1;-1;1;1;1;1;0;1;-1;-1;1;
    1;-1;1;-1;1;-1;-1;-1;-1;-1;1;1;-1;-1;1;-1;
    1;-1;1;1;1;1];
LTFpriorIFFT = [zeros(6,1);seqLong;zeros(5,1)];
LTF_ifft = ifft(ifftshift(LTFpriorIFFT,1))*(1/sqrt(52))*64;
LTFifft_double = [LTF_ifft;LTF_ifft];
LTFifft_full = [LTFifft_double(97:128);LTFifft_double];

%% Cross-correlation
RX_Data_picked  = dataset_wifi(sampleNum * (nnn-1) + 1: sampleNum * (nnn) );
% RX_Data_picked  = RX_Data((nnn-1) * sampleNum +1:nnn * sampleNum );
RX_Data_picked = rx;
[r,~] = xcorr(RX_Data_picked,LTFifft_full);
cross_correlation = r(length(RX_Data_picked):end-length(LTFifft_full)+1);
squared_absolute = abs(cross_correlation).^2;
figure(2)
plot(squared_absolute);

[~,max_index] = max(squared_absolute);
if max_index>160
    offset_zzf = max_index - 160;
else
    offset_zzf = max_index -160 + sampleNum;
end

 offsetData = RX_Data ( offset_zzf + sampleNum * (nnn-1) : sampleNum * nnn+  offset_zzf -1 );
 % 30.72
 frame_sync = offsetData(1: 2960);

%% WLAN Function Synchronization
% % NOTE 1 : The output metric is the squared absolute value rather than
% %just absolute value.
% % NOTE 2: Offset +1 is the start of the frame rather than offset, offset
% % just means the deviated number of samples
[startOffset,meric] = wlanSymbolTimingEstimate(RX_Data_picked,'CBW20'); %% rxSig  RX_Data WIFI_Waveform
figure(1)
plot(meric)

isequal(startOffset, offset_zzf -1)





