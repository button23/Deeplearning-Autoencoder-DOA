function frame_sync = IPS_synchronization(input_wlan, cfgVHT,  totalNum)
%% RF data input format conversion
input_reshape = reshape(input_wlan,2,numel(input_wlan)/2);

% Imaginary and Real part separation
imag_input = input_reshape(1,:).';
real_input = input_reshape(2,:).';

RX_Data = complex(real_input, imag_input);

if cfgVHT.MCS ==8
    sampleNum = 2960;
else
    sampleNum = 960;
end

%% WLAN Function Synchronization
% NOTE 1 : The output metric is the squared absolute value rather than
%just absolute value.
% NOTE 2: Offset +1 is the start of the frame rather than offset, offset
% just means the deviated number of samples
% [startOffset,meric] = wlanSymbolTimingEstimate(RX_Data,'CBW20'); %% rxSig  RX_Data WIFI_Waveform
% figure(1)
% plot(meric)


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
RX_Data_picked  = RX_Data(1: sampleNum );
% RX_Data_picked  = RX_Data((nnn-1) * sampleNum +1:nnn * sampleNum );

[r,~] = xcorr(RX_Data_picked,LTFifft_full);
cross_correlation = r(length(RX_Data_picked):end-length(LTFifft_full)+1);
squared_absolute = abs(cross_correlation).^2;
% figure(2)
% plot(squared_absolute);

[~,max_index] = max(squared_absolute);
if max_index>160
    offset_zzf = max_index - 160;
else
    offset_zzf = max_index -160 + sampleNum;
end

 offsetData = RX_Data ( offset_zzf : totalNum * sampleNum + offset_zzf - 1);
 frame_sync = reshape(offsetData, sampleNum, totalNum);
end
