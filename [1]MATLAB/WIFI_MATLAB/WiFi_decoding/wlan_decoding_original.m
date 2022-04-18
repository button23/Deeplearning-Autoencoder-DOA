%% Initialization
% clear;
% clc;
% close all;
self_scaling = 2^13;
BM_scaling = 2^8;
%% Data Convertion
% VarName1_reshape = reshape(VarName1,2,length(VarName1)/2);
% imag_data = VarName1_reshape(2,:);
% real_data = VarName1_reshape(1,:);
% waveform_data = real_data + 1j*imag_data;
rxSig = frame_sync;  %%frame_sync  synchFrame

%% wlan waveform
cfgVHT = wlanVHTConfig('ChannelBandwidth','CBW20','MCS',8);
% cfgVHT
% Coded bits per single carrier for each spatial steam
Nbpscs = 8;
R = 3/4;
% Coded bits per symbol (52 subcarrier * 8bits(256 QAM))
Ncbps = 52 * Nbpscs;
% Data bits per OFDM symbol (416 * 3/4 = 312 bits)
Ndbps = Ncbps * R;
% Used for scaling Matlab data in order to move to DSP
scalingFactor = 2^11;
% Default value is 1024
cfgVHT.APEPLength = 1024;
% Create random bit stream
rng(25)
psdu = randi([0 1],cfgVHT.PSDULength*8,1);

% Convert psdu into ppdu
wholeWaveform = wlanWaveformGenerator(psdu,cfgVHT);
wholeWaveform_scaling = scalingFactor* wholeWaveform;

ind_Data = wlanFieldIndices(cfgVHT,'VHT-Data');
ind_VHTLTF = wlanFieldIndices(cfgVHT,'VHT-LTF');

ind_STF = wlanFieldIndices(cfgVHT,'L-STF');
ind_LTF = wlanFieldIndices(cfgVHT,'L-LTF');
ind_LSIG = wlanFieldIndices(cfgVHT,'L-SIG');

data_part_transmit = wholeWaveform(ind_Data(1):ind_Data(2),:);

% psdu = int8(psdu);
% y = wlanVHTData(psdu,cfgVHT);  % 93 default     2160= 801:2960
[sigB,bits_sigB] = wlanVHTSIGB(cfgVHT);

%% CRC of SIGB
sigB_20bit = bits_sigB(1:20);
window = ones(8,1);
for n=1:length(sigB_20bit)
    a = xor(window(1),sigB_20bit(n));
    window(1:5) = window(2:6);
    window(6) = xor(window(7),a);
    window(7) = xor(window(end),a);
    window(end) = a;
end
crc = ~window;

%% Add Service Field and tail bits
VHTData_Service = zeros(16,1);
VHTData_Service(9:16) = crc;
% VHTData_Service(9:16) = [1;1;1;1;1;1;1;1];

% Calculate the number of OFDM symbols
Nsym = ceil((cfgVHT.PSDULength*8+16+6)/Ndbps);
% Calculate the padding bits
Npad = Nsym*Ndbps - cfgVHT.PSDULength*8 - 16 - 6;
% Insert padding bits and service field into data field
VHTData = [VHTData_Service;psdu;zeros(Npad,1);zeros(6,1)];

% VHTData = int8(VHTData);

%% Scrambler
% scramInit = 93;
% scrambleData = wlanScramble(VHTData,scramInit);
%%Replace the srambled tail bits with zero bits.
% scrambleData =[1;0;1;1;1;0;1;scrambleData_(8:end-6);zeros(6,1)];

% Initial value 93
scrambler =[ 1 0 1 1 1 0 1];
scrambleData = zeros(length(VHTData),1);
for scr = 1:length(VHTData)
    tem = xor(scrambler(1),scrambler(4));
    scrambler(1:end-1) = scrambler(2:end);
    scrambler(end) = tem;
    scrambleData(scr) = xor(tem,VHTData(scr));
end
scrambleData =[scrambleData(1:end-6);zeros(6,1)];

%% BCC encoder
BCCencodedData = wlanBCCEncode(scrambleData,R);

%% BCC Interleaver
kkk = (0:Ncbps-1)';
interleaver_data = zeros(Ncbps,Nsym);
iii = zeros(Ncbps,1);
jjj = zeros(Ncbps,1);
s = max(1,Nbpscs/2);

% cbw = 'CBW20';
% interleaver_data_test = wlanBCCInterleave(BCCencodedData,'VHT',Ncbps,cbw);

permutation_first = zeros(Ncbps,1);
permutation_second = zeros(Ncbps,1);


for internn = 1:Nsym
    first_Ncbps = BCCencodedData(1+Ncbps*(internn-1):Ncbps*internn);
    
    % first index permutation
    iii = Ncbps/13 *mod(kkk,13)+floor(kkk/13);
    
    for cc = 1:length(kkk)
        permutation_first(iii(cc)+1) = first_Ncbps(cc);
    end
    % second index  permutation
    jjj = s*floor(kkk/s)+mod((kkk+Ncbps-floor(kkk*13/Ncbps)),s);
    
    for cc = 1:length(kkk)
        permutation_second(jjj(cc)+1) = permutation_first(cc);
    end
    interleaver_data(:,internn) = permutation_second;
end
% interleaver_data = int8(interleaver_data);

%% Constellation Mapper
mappedData = wlanConstellationMap(interleaver_data,Nbpscs);
% modTx = qammod(interleaver_data,64,'InputType','bit');

%% polarity sequence
% Used for pilot bit (served as sign)
pv=[1;1;1;1; -1;-1;-1;1; -1;-1;-1;-1; 1;1;-1;1; -1;-1;1;1; -1;1;1;-1; 1;1;1;1; 1;1;-1;1;
    1;1;-1;1; 1;-1;-1;1; 1;1;-1;1; -1;-1;-1;1; -1;1;-1;-1; 1;-1;-1;1; 1;1;1;1; -1;-1;1;1;
    -1;-1;1;-1; 1;-1;1;1; -1;-1;-1;1; 1;-1;-1;-1; -1;1;-1;-1; 1;-1;1;1; 1;1;-1;1; -1;1;-1;1;
    -1;-1;-1;-1; -1;1;-1;1; 1;-1;1;-1; 1;1;1;-1; -1;1;-1;-1; -1;1;1;1; -1;-1;-1;-1; -1;-1;-1];

Data_pilot = zeros(64,Nsym);
pil=zeros(4,4);
% The following four sequences are repeatedly used for pilot bit of OFDM
% symbol (the sequence should be multiplied by the corresponding polarity bit)
pil(:,1) = [1,1,1,-1];
pil(:,2) = [1,1,-1,1];
pil(:,3) = [1,-1,1,1];
pil(:,4) = [-1,1,1,1];

%% insert pilots
for nnss = 1:Nsym
    %     data_prepilot = mappedData(1+52*(nnss-1):52*nnss);
    data_prepilot = mappedData(:,nnss);
    pilo = pv(nnss+4)*pil(:,mod((nnss-1),4)+1);
    Data_pilot(:,nnss) = [zeros(4,1);data_prepilot(1:7);pilo(1);data_prepilot(8:20);pilo(2);data_prepilot(21:26);0;data_prepilot(27:32);pilo(3);data_prepilot(33:45);pilo(4);data_prepilot(46:52);zeros(3,1)];
end

%% IFFT
% k = -32:1:31;
% k =  k' ;
% vhtdata_ifft = zeros(64,Nsym);
% for ifftnn = 1:Nsym
%     for kk = 1:64
%             yy =  Data_pilot(:,ifftnn) .* exp(j*2*pi*k*(kk-1)/64);
%             vhtdata_ifft(kk,ifftnn) = 1 / sqrt(56) * sum(yy) ;
%     end
% end

% If using ifft function, we should ifftshift first and then ifft.
% 1/sqrt(56) is the scaling for WIFI, and 64 is because usually in fft,
% there would be 64 scaling, but here is the opposite, in ifft there would
% be scaling 64, therefore, we need to multiply 64 in ifft
vhtdata_ifft =(1/sqrt(56))* 64 * ifft(ifftshift(Data_pilot,1));

%% Add CP
cp_data = zeros(80,Nsym);
for cpnn = 1:Nsym
    cp_data(:,cpnn) = [ vhtdata_ifft(49:64,cpnn);vhtdata_ifft(:,cpnn)];
end
transmittt = reshape(cp_data,[80*Nsym,1]);

% % To verity if the self-made data field is the same as the RMC waveform.
% It turned out that only the first sample of each OFDM symbol is
% different due to the pulse shape.
% figure
% plot(abs(data_part_transmit))
% hold on
% plot(abs(transmittt))
% legend("Matlab","Self-made")
% xlabel('Data Field Waveform')
% ylabel('Amplitute')

%% Create a phase and frequency offset object and introduce a 2 Hz frequency offset.
tgacChan = wlanTGacChannel('SampleRate',20e6,'ChannelBandwidth','CBW20','DelayProfile','Model-C');

chNoise = comm.AWGNChannel('NoiseMethod','Signal to noise ratio (SNR)',...
    'SNR',35);
pfOffset = comm.PhaseFrequencyOffset('SampleRate',20e6,'FrequencyOffsetSource','Input port');
% rxSig = chNoise(tgacChan(wholeWaveform));  %% wholeWaveform_scaling  wholeWaveform
% rxSig = chNoise(wholeWaveform_scaling);
% rxSig = pfOffset(rxSig,10000);
% rxSig = tgacChan(chNoise(pfOffset(wholeWaveform)));
% rxSig = wholeWaveform;

%% Data recover (WLAN function)
% recBits = wlanVHTDataRecover(transmittt,ones(56,1),0,cfgVHT);
% numErrs = biterr(psdu,recBits)

%% Coarse CFO Estimation
% Using the last five identical parts to do coarse frequency offset
% estimation
stf_part = rxSig(ind_STF(1):ind_STF(2),:);
fisrtPartStf = stf_part(81:144);
secondPartStf = stf_part(97:160);
% Caluculate the phase shift between every sample
phase_STF = 1/16 * angle(fisrtPartStf'*secondPartStf);
% the calculated phase is in the unit of radian,change it to frequency.
Ts = 1/(2*10^7);
CFO_STF = phase_STF/(2*pi*Ts)

%% Coarse CFO Correlation
% Compensate only the LTF, and the index is set to be zero instead of
% the initial index in STF to symplify the process
% because LTF would do later calibration any way
mmC = 0 : 127;
% Be caucious, this line can not be run multiple times, because the rxSig
% would be overwritten
% rxSig(ind_LTF(1)+32:ind_LTF(2)) = rxSig(ind_LTF(1)+32:ind_LTF(2)) .* exp(-1i * mmC' * phase_STF);


%% Fine CFO Estimation
Ts = 1/(2*10^7);
ltf_part = rxSig(ind_LTF(1):ind_LTF(2),:);
%
first_LTF = ltf_part(33:96);
% % For dsp data
% real_ltf_first = int16(real(first_LTF)*1000);
% imag_ltf_first = int16(imag(first_LTF)*1000);
%
second_LTF = ltf_part(97:160);
% % For dsp data
% real_ltf_second = int16(real(second_LTF)*1000);
% imag_ltf_second = int16(imag(second_LTF)*1000);
%
phase_LTF = 1/64 * angle(first_LTF'*second_LTF);
CFO_LTF = phase_LTF/(2 * pi * Ts)
% CFO_LTF = 130;
% phase_LTF = CFO_LTF *(2 * pi * Ts);
% % Total carrier frequency offset phase combined with coarse CFO estimation and
% % Fine CFO estimation
total_phase = phase_LTF + phase_STF
% total_cfo = CFO_LTF + CFO_STF;

%% Wlan Funciton to estimate CFO
% WIFI function
% foffset1 = wlanCoarseCFOEstimate(stf_part,'CBW20')
% rxSig= pfOffset(rxSig,-foffset1);
% ltf_part = rxSig(ind_LTF(1):ind_LTF(2),:);
% foffset2 = wlanFineCFOEstimate(ltf_part,'CBW20')
% rxSig= pfOffset(rxSig,-foffset2);

%% Fine CFO Correlation
% mmF = 128:128+ind_Data(2)-ind_LTF(2)-1;
mmF = 128:128+ind_Data(2)-ind_LTF(2)-1;

mmF = double(mmF');

% If do Fine Compensation at the beginning, the error rate is less and the
% constallation map is  neater.
% rxSig(ind_LTF(2)+1:end) = rxSig(ind_LTF(2)+1:end) .* exp(-j * mmF * total_phase);
% rxSig(ind_LTF(2)+1:end) = rxSig(ind_LTF(2)+1:end) .* exp(-j * mmF * phase_LTF);

%% DSP RF Data Preparation
rf_data_int16 = int16(rxSig);
rf_data_int16_real = real(rf_data_int16);
rf_data_int16_imag = imag(rf_data_int16);
rf_data_ = [rf_data_int16_imag';rf_data_int16_real'];
rf_data_combine = reshape(rf_data_,2960*2,1);

rf_data_combine_4frames = repmat(rf_data_combine,[4,1]);

% add offset to this combined four frame
rf_data_combine_4frames_offset = [rf_data_combine_4frames(1560:end);rf_data_combine_4frames(1:1559)];

%% Remove CP (L-LTF)
% cp_remove_llf = ltf_part(33:end);
% vhtllf_rem_cp = reshape(cp_remove_llf,[64,2]);

%% FFT( L-LTF )
% k2 = 1:1:64;
% k2 = k2';
% l_llf_fft = zeros(64,2);
% for fftnn = 1:2
%     for kk = 1:64
%             yy =  sqrt(52) *vhtllf_rem_cp(:,fftnn) .* exp(-j*2*pi*(kk-33)*(k2-1)/64);
%             l_llf_fft(kk,fftnn) = 1/64 *  sum(yy) ;
%     end
% end


%% Remove zero padding ( L-LTF )
% l_llf_af_pilots = [l_llf_fft(7:32,:);l_llf_fft(34:59,:)];

%% Channel Estimation (L-LTF)
% seqLong_L = [1;1;-1;-1;1;1;-1;1;-1;1;1;1;1;1;1;
%             -1;-1;1;1;-1;1;-1;1;1;1;1;1;-1;-1;1;
%             1;-1;1;-1;1;-1;-1;-1;-1;-1;1;1;-1;-1;1;-1;
%             1;-1;1;1;1;1];
% H_L = (l_llf_af_pilots(:,1)+l_llf_af_pilots(:,2))/2 .* seqLong_L;

%% Remove CP(VHT-LTF)
vhtltf_part = rxSig(ind_VHTLTF(1):ind_VHTLTF(2),:);
cp_remove_vhtltf = vhtltf_part(17:end);

% For DSP input
% fft_vhtltf_real = int16(real(cp_remove_vhtltf*sqrt(56)/64));
% fft_vhtltf_imaginary = int16(imag(cp_remove_vhtltf*sqrt(56)/64));

% Without scaling sqrt(56)/64 (when doing channel compensation, because both denominator and
%numerator contain this term; therefore, it will be cancelled out)
% fft_vhtltf_real = int16(real(cp_remove_vhtltf*scalingFactor));
% fft_vhtltf_imaginary = int16(imag(cp_remove_vhtltf*scalingFactor));
% fft_vhtltf_combine = [fft_vhtltf_real;fft_vhtltf_imaginary];

%% FFT( VHT-LTF )
k2 = 1:1:64;
k2 = k2';
vht_llf_fft = zeros(64,1);

for kk = 1:64
    yy =  sqrt(56) *cp_remove_vhtltf(:) .* exp(-j*2*pi*(kk-33)*(k2-1)/64);
    vht_llf_fft(kk) = 1/64 *  sum(yy) ;
end

% dsp_test_vhtltf = fft(cp_remove_vhtltf*sqrt(56)/64);

% Without scaling sqrt(56)/64 (when doing channel compensation, because both denominator and
%numerator contain this term; therefore, it will be cancelled out)
% set the static scaling of FFTC to be 1/4, so that even without
% sqrt(56)/64, the data won't exceed 32768

% Matlab Data (needs scaling)
% dsp_test_vhtltf = fft(cp_remove_vhtltf*scalingFactor);
% dsp_test_vhtltf_fftshift = double(fftshift(dsp_test_vhtltf,1));

% RF data (with scaling)
dsp_test_vhtltf = fft(cp_remove_vhtltf);
dsp_test_vhtltf_fftshift = double(fftshift(dsp_test_vhtltf,1));

% To check the dynamic scaling of FFTC in DSP
dsp_test_vhtltf_fftshift_1 =  int32(dsp_test_vhtltf_fftshift ) ;
dsp_test_vhtltf_fftshift_2 =  int32(dsp_test_vhtltf_fftshift * 2) ;
dsp_test_vhtltf_fftshift_1_2 =  int32(dsp_test_vhtltf_fftshift / 2) ;
dsp_test_vhtltf_fftshift_1_4 =  int32(dsp_test_vhtltf_fftshift / 4) ;
dsp_test_vhtltf_fftshift_1_8 =  int32(dsp_test_vhtltf_fftshift / 8) ;

%% Remove zero padding ( VHT-LTF )
% This one is without scalingFactor(original data)
vht_llf_af_pilots = [vht_llf_fft(5:32);vht_llf_fft(34:61)];
% This one is with scalingFactor(Fpr DSP data)
vht_llf_af_pilots_scaling = [dsp_test_vhtltf_fftshift(5:32);dsp_test_vhtltf_fftshift(34:61)];

%% Channel Estimation ( VHT-LTF )
seqLong_VHT = [1;1;1;1;-1;-1;1;1;-1;1;-1;1;1;1;1;1;1;
    -1;-1;1;1;-1;1;-1;1;1;1;1;1;-1;-1;1;
    1;-1;1;-1;1;-1;-1;-1;-1;-1;1;1;-1;-1;1;-1;
    1;-1;1;1;1;1;-1;-1];
H_VHT = vht_llf_af_pilots .* seqLong_VHT;
% figure
% plot(abs(H_VHT))

% Channel estimation with scaling factor.(DSP data)
H_VHT_scaling = vht_llf_af_pilots_scaling .*seqLong_VHT;

%% Remove CP(data)
data_part = rxSig(ind_Data(1):ind_Data(2),:);
res_data_bf_CP = reshape(data_part,[80,Nsym]);
data_cp_remove= res_data_bf_CP(17:80,:);


%% FFT (data)
% For DSP data input ( DSP Purpose )
% fft_data_real = int16(real(scalingFactor*data_cp_remove*sqrt(56)/64));
% fft_data_imaginary = int16(imag(scalingFactor*data_cp_remove*sqrt(56)/64));

% Without scaling sqrt(56)/64 (when doing channel compensation, because both denominator and
%numerator contain this term; therefore, it will be cancelled out)
fft_data_real = int16(real(data_cp_remove*scalingFactor));
fft_data_imaginary = int16(imag(data_cp_remove*scalingFactor));
dsp_fft_data_combine = [fft_data_real;fft_data_imaginary];
% Put the following one directly into fft_testData.h
dsp_fft_data_combine_one_column = reshape(dsp_fft_data_combine,128*Nsym,1);

k2 = 1:1:64;
k2 = k2';
vhtdata_fft = zeros(64,Nsym);
for fftnn = 1:Nsym
    for kk = 1:64
        yy =  sqrt(56) *data_cp_remove(:,fftnn) .* exp(-j*2*pi*(kk-33)*(k2-1)/64);
        vhtdata_fft(kk,fftnn) = 1/64 *  sum(yy) ;
    end
end

% First fft and then fftshift, the sequence
% is reversed compared to ifft on the transmitter side
% dsp_test_data = fft(data_cp_remove*sqrt(56)/64);

% Without scaling sqrt(56)/64 (when doing channel compensation, because both denominator and
%numerator contain this term; therefore, it will be cancelled out)
% set the static scaling of FFTC to be 1/4, so that even without
% sqrt(56)/64, the data won't exceed 32768

% Matlab Data (needs scaling)
% dsp_test_data = fft(data_cp_remove*scalingFactor);
% dsp_test_fftshift = fftshift(dsp_test_data,1);

% RF data (with scaling)
dsp_test_data = fft(data_cp_remove);
dsp_test_fftshift = fftshift(dsp_test_data,1);


% To verify the result of FFTC with dynamic scaling
dsp_test_fftshift_1 = int32(dsp_test_fftshift);
dsp_result_fft_scale_2 = int32(dsp_test_fftshift*2);
dsp_result_fft_scale_4 = int32(dsp_test_fftshift*4);
dsp_result_fft_scale_8 = int32(dsp_test_fftshift*8);

dsp_result_fft_scale_1_2 =  int32(dsp_test_fftshift / 2) ;
dsp_result_fft_scale_1_4 =  int32(dsp_test_fftshift / 4) ;
dsp_result_fft_scale_1_8 =  int32(dsp_test_fftshift / 8) ;

% plot to check whether the scaling factor is appropriate (if any real part of imaginary part
% exceeds the limititon of int16 -32768~32765)
% figure
% plot(abs(real(dsp_test_fftshift)))
% hold on
% plot(abs(imag(dsp_test_fftshift)))
% legend("Real Part","Imaginary Part")
% xlabel('Data Field after FFTC ')
% ylabel('Amplitute')



%% Remove zero padding
% zeropadding_remove_origin = [vhtdata_fft(5:32,:);vhtdata_fft(34:61,:)];
zeropadding_remove = [dsp_test_fftshift(5:32,:);dsp_test_fftshift(34:61,:)]; % With scaling

%% Residual Frequency Offset Estimation
% % pilot_data_origin = [zeropadding_remove_origin(8,:);zeropadding_remove_origin(22,:);zeropadding_remove_origin(35,:);zeropadding_remove_origin(49,:)];
% pilot_data = [zeropadding_remove(8,:);zeropadding_remove(22,:);zeropadding_remove(35,:);zeropadding_remove(49,:)]; % With scaling
% H_pilot = [H_VHT(8);H_VHT(22);H_VHT(35);H_VHT(49)];
% H_pilot_scaling = [H_VHT_scaling(8);H_VHT_scaling(22);H_VHT_scaling(35);H_VHT_scaling(49)];
% 
% residual_offset = zeros(Nsym,1);
% residual_offset_amplitude = zeros(Nsym,1);
% 
% for resnn =1:Nsym
%     pilo_ = pv(resnn+4)*pil(:,mod((resnn-1),4)+1);
%     residual_offset(resnn) = angle(pilot_data(:,resnn)'*(H_pilot.*pilo_));
% end
% 
% 
% %Test
% compensated_real_residual = zeros(56,Nsym);
% compensated_imag_residual = zeros(56,Nsym);
% for resnn =1:Nsym
%     pilo_ = pv(resnn+4)*pil(:,mod((resnn-1),4)+1);
%     residual_offset_amplitude(resnn) = pilot_data(:,resnn)'*(H_pilot.*pilo_)/4;
% end
% residual_real = real(residual_offset_amplitude);
% residual_imag = imag(residual_offset_amplitude);
% 
% residual_amp_ = abs(residual_offset_amplitude);
% 
% 
% nom_residual = (residual_real .* residual_real + residual_imag.* residual_imag)/ 2^14;
% 
% for estnn = 1:Nsym
%     aaa = real(zeropadding_remove(:,estnn)) .* residual_real(estnn);  %zeropadding_remove has already been factored, /4 is dynamic FFTC scaling
%     bbb = imag(zeropadding_remove(:,estnn)) .* residual_imag(estnn);
%     ccc = imag(zeropadding_remove(:,estnn)) .* residual_real(estnn);
%     ddd = real(zeropadding_remove(:,estnn)) .* residual_imag(estnn);
%     
%     %     compensated_real(:,estnn) = ((aaa + bbb)./ nom);
%     %     compensated_imag(:,estnn) = ((ccc - ddd)./ nom);
%     compensated_real_residual(:,estnn) = ((aaa + bbb))/ 2^14;
%     compensated_imag_residual(:,estnn) = ((ccc - ddd))/ 2^14;
% end
% 
% compensated_full_residual = compensated_real_residual + i*compensated_imag_residual;
% for estnn = 1:Nsym
%     need_to_divide = compensated_full_residual(:,estnn)./nom_residual(estnn);
% end


%% Channel Compensation and Residual Frequency Offset compensation
aft_estcom = zeros(56,Nsym);
% The following way can not be used for DSP, because of the scaling of the
% result. Numerator and denominator have the same degree which results in
% the 0 degree result
for estnn = 1:Nsym
    %     aft_estcom(:,estnn) =  zeropadding_remove(:,estnn)./(H_VHT) * exp(1i*residual_offset(estnn));
    %     aft_estcom(:,estnn) =  (zeropadding_remove(:,estnn))./vht_llf_scaling_est ;   %%%    H_VHT is not even scaled.  vht_llf_scaling_est H_VHT
    %     aft_estcom(:,estnn) =  (zeropadding_remove(:,estnn))./(H_VHT) ;   %%  This one can be used to observe the CFO
        aft_estcom(:,estnn) =  (zeropadding_remove(:,estnn))./(H_VHT_scaling) ;   %%  This one can be used to observe the CFO
%     aft_estcom(:,estnn) =  (zeropadding_remove_origin(:,estnn))./(H_VHT) ;   %%  This one can be used to observe the CFO
    
%         aft_estcom(:,estnn) =  zeropadding_remove(:,estnn)./(H_VHT_scaling) * exp(1i*residual_offset(estnn));
    %     aft_estcom(:,estnn) =  (zeropadding_remove(:,estnn));
end



% Test new way for dsp (Instead of dividing the denominator, multiply it
% with sqrt(1/170), which is the coefficient in 256QAM de-modulation )
compensated_real = zeros(56,Nsym);
compensated_imag = zeros(56,Nsym);
test_vht_real = real(H_VHT_scaling ) ;
test_vht_imag = imag(H_VHT_scaling ) ;
% The reason why 2^14 is chosen is out of experiment( appropriate scaling)
% divided by 16 is due to divide the square of fftc dynamic scaling 8,
% 8*8=64
nom = (test_vht_real .* test_vht_real + test_vht_imag.* test_vht_imag) / self_scaling;
nom_max = int16(max(nom));


for estnn = 1:Nsym
    aaa = real(zeropadding_remove(:,estnn)) .* test_vht_real;
    bbb = imag(zeropadding_remove(:,estnn)) .* test_vht_imag;
    ccc = imag(zeropadding_remove(:,estnn)) .* test_vht_real; % bc
    ddd = real(zeropadding_remove(:,estnn)) .* test_vht_imag; % ad
    compensated_real(:,estnn) = ((aaa + bbb)) / self_scaling;
    compensated_imag(:,estnn) = ((ccc - ddd)) / self_scaling;
    %     compensated_real(:,estnn) = ((aaa + bbb)) / (8*8) / 2^14* exp(1i*residual_offset(estnn));
    %     compensated_imag(:,estnn) = ((ccc - ddd)) / (8*8) / 2^14* exp(1i*residual_offset(estnn));
end
compensated_real_max = int16(max(compensated_real));

%% Remove pilots
% pilot_remove =[zeropadding_remove(1:7,:);zeropadding_remove(9:21,:);zeropadding_remove(23:34,:);zeropadding_remove(36:48,:);zeropadding_remove(50:56,:)];
% Following operaiton without practial usage, without scaling, not suited
% for DSP
pilot_remove =[aft_estcom(1:7,:);aft_estcom(9:21,:);aft_estcom(23:34,:);aft_estcom(36:48,:);aft_estcom(50:56,:)];


% DSP test 
pilot_remove_real_test =[compensated_real(1:7,:);compensated_real(9:21,:);compensated_real(23:34,:);compensated_real(36:48,:);compensated_real(50:56,:)];
pilot_remove_imag_test =[compensated_imag(1:7,:);compensated_imag(9:21,:);compensated_imag(23:34,:);compensated_imag(36:48,:);compensated_imag(50:56,:)];
pilot_remove_nom = [nom(1:7,:);nom(9:21,:);nom(23:34,:);nom(36:48,:);nom(50:56,:)];
pilot_remove_nom_1_2 = int16(pilot_remove_nom);
A_scaling  =  sqrt(1/170) * pilot_remove_nom_1_2;

% The following one is without any scaling, it goes back to the origin.
% This result should be the same as pilot_remove
pilot_remove_test = (pilot_remove_real_test+i*pilot_remove_imag_test)./pilot_remove_nom ;
% The following is the real data case (with scaling)
pilot_remove_test_scaling = int16((pilot_remove_real_test+i*pilot_remove_imag_test)) ;
pilot_remove_nom_scaling = pilot_remove_nom./2^5;
pilot_remove_test_scaling_withConMaP = (pilot_remove_real_test+i*pilot_remove_imag_test)./pilot_remove_nom_scaling ;
pilot_remove_test_scaling_one = reshape(pilot_remove_test_scaling_withConMaP,52*27,1);
% Residual CFO compensation
% pilot_remove_test_with_residule_com = zeros(52,Nsym);
% for estnn = 1:Nsym
% pilot_remove_test_with_residule_com(:,estnn) = pilot_remove_test(:,1)* exp(1i*residual_offset(estnn));
% end
%
% figure(5)
% pilot_remove_test_ =(reshape(pilot_remove_test,[52*Nsym,1]));
% scatter(real(pilot_remove_test_),imag(pilot_remove_test_))

% plot the result of compensated data
figure(6)
xlim([-1.5 1.5])
ylim([-1.5 1.5])
% scatterplot(reshape(pilot_remove,[52*Nsym,1]))
pilot_remove_ =(reshape(pilot_remove,[52*Nsym,1]));
% scatter(real(pilot_remove_),imag(pilot_remove_),'filled')
% scatter(real(pilot_remove_),imag(pilot_remove_),28,'MarkerEdgeColor',[0 .1 .5],...
%               'MarkerFaceColor',[0 .5 .7],...
%               'LineWidth',1.0)
sz = 25;
scatter(real(pilot_remove_),imag(pilot_remove_),sz,'*')
% scatter(real(pilot_remove_test_scaling_one),imag(pilot_remove_test_scaling_one),sz,'*')
% scatter(real(RX_Data),imag(RX_Data),sz,'*')


xlabel('In-phase')
ylabel('Quadrature')
title('Received 256 QAM Constellation')

%% Constellation Demapper(algorithm)
llrAggregate = zeros(Ncbps, Nsym);
llr = zeros(8,1);
%if there is no compensation(or during compensation Scaling factor is not divided, then here should be added with scaling factor)
%  A = sqrt(1/170)*scalingFactor;
A = sqrt(1/170);

qamBits= zeros(Nbpscs,52);
for oo = 1  :  Nsym
    for ll = 1: 52
%         A = A_scaling(ll);
%         rr = real(pilot_remove_test_scaling(ll,oo));
%         im = imag(pilot_remove_test_scaling(ll,oo));
                    rr = real(pilot_remove(ll,oo));
                    im = imag(pilot_remove(ll,oo));
        llr(1) = rr ;   %% no need for /4 as long as zeropadding_remove does not divide anything
        llr(2) = - abs(llr(1)) + 8 * A;
        llr(3) = - abs( llr(2) ) + 4 * A;
        llr(4) = - abs( llr(3)) + 2 * A;
        
        llr(5) = im ;
        llr(6) = - abs(llr(5)) + 8 * A;
        llr(7) = - abs( llr(6)) + 4 * A;
        llr(8) = - abs( llr(7)) + 2 * A;
        
        %             qamBits(:,ll) = llr * (abs(H_VHT(ll)))^2 ;
        qamBits(:,ll) = -llr ;
    end
    % Minus sign to inverse the LLR value
    llrAggregate(:,oo) = reshape(qamBits,Ncbps,1);
end
llrAggregate_scale = int16(llrAggregate);
demod_Algorithm = llrAggregate < 0;
% numError = biterr(demod_Algorithm,interleaver_data)
% ber = numError / (Ncbps*Nsym)

%% Constellation Demapper
% demappedDataRef = wlanConstellationDemap(pilot_remove,0,Nbpscs,'soft');

%% BCC Deinterleaver
% cbw = 'CBW20';
% demappedDataRef_one_column = reshape(demappedDataRef,Ncbps*Nsym,1);
% deinterleaver_data_refer = wlanBCCDeinterleave(demappedDataRef_one_column,'VHT',Ncbps,cbw);

%% BCC Deinterleaver (code)
kkk_de = (0:Ncbps-1)';
deinterleaver_data = zeros(Ncbps,Nsym);
iii_de = zeros(Ncbps,1);
jjj_de = zeros(Ncbps,1);
s = max(1,Nbpscs/2);

permutation_first = zeros(Ncbps,1);
permutation_second = zeros(Ncbps,1);

for internn = 1:Nsym
    first_Ncbps = llrAggregate(:,internn);
    
    % first index permutation]
    iii_de = s*floor(kkk_de/s)+mod((kkk_de+floor(kkk_de*13/Ncbps)),s);
    
    for cc = 1:length(kkk_de)
        permutation_first(iii_de(cc)+1) = first_Ncbps(cc);
    end
    
    % second index  permutation
    %     jjj_de = Ncbps/13 *mod(kkk_de,13)+floor(13*kkk_de/);
    jjj_de = 13*kkk_de - (Ncbps-1)*floor(kkk_de*13/Ncbps);
    
    for cc = 1:length(kkk_de)
        permutation_second(jjj_de(cc)+1) = permutation_first(cc);
    end
    deinterleaver_data(:,internn) = permutation_second;
end
deinter_data_reshape = reshape(deinterleaver_data,Ncbps*Nsym,1);


%% De-Puncturing
addzeros = zeros(2,1);
after_puncturing = zeros(length(VHTData)*2,1);
numAddZeros = Ncbps*Nsym /4-1;
after_puncturing(1:5) = [deinter_data_reshape(1:3);addzeros];
after_puncturing(end) = deinter_data_reshape(end);
for naz = 1:numAddZeros
    index_naz_des = 6 * naz;
    index_naz_orig = 4 * naz;
    
    after_puncturing(index_naz_des:index_naz_des+5) = [deinter_data_reshape(index_naz_orig:index_naz_orig+3);addzeros];
    [m,n] = size(after_puncturing);
    puncturing_reshape = reshape(after_puncturing,[(m*n), 1]);
end

%% Calculate Branch Metrics (BM)
BM = zeros(length(puncturing_reshape),1);
for ii=1:2:length(puncturing_reshape)
    BM(ii) = puncturing_reshape(ii)+puncturing_reshape(ii+1);
    BM(ii+1) = puncturing_reshape(ii)-puncturing_reshape(ii+1);
end

% [maxi,ind] = max(abs(BM));
% if BM(ind) < 0
%     BM_scale = BM / maxi * 64;
% else
%     BM_scale = BM / maxi * 63;
% end
BM_scale = BM /BM_scaling;
[maxi,ind] = max(abs(BM_scale));

BM_scale_int8 = int8(BM_scale);

BM_scale_reshape = reshape(BM_scale_int8,624,Nsym);

%% BCC Decoder
% decodedData_reference = wlanBCCDecode(deinter_data_reshape,R,'soft');

% decodedData_stuff_zero = [decodedData_reference;zeros(8,1)]; %% Stuff 0 to make the number be the multiple of 32
%
% decodedData_reference_trans = reshape(decodedData_stuff_zero,32,88); %%(2808+8)/32 = 88
% decodedData_reference_trans_ = decodedData_reference_trans';
% decodedData_deci = bi2de(decodedData_reference_trans_);

% WIFI decoder function
decodedData = wlanBCCDecode(puncturing_reshape,1/2,'soft');

% Viterbi decoder setting
trellis = poly2trellis(7,[133 171]);
tbl = 32;
rate = 1/2;
% Viterbi decode the demodulated data
% dataHard = vitdec(rxDataHard,trellis,tbl,'cont','hard');
dataSoft = vitdec(puncturing_reshape,trellis,tbl,'trunc','unquant');

% below is for DSP result confirmation (decisions of VCP2)
dataSoft_stuff_zero = [dataSoft;zeros(24,1)];
dataSoft_re = reshape(dataSoft_stuff_zero,32,numel(dataSoft_stuff_zero)/32);
dataSoft_re_ = dataSoft_re';
dataSoft_deci=int64(bi2de(dataSoft_re_));

%% Descrambler
% scramInit = 93;
% receive_descramble = wlanScramble(decodedData,scramInit);
% errorDescram = sum(abs(int8(receive_descramble)-int8(VHTData)))

scrambler =[ 1 0 1 1 1 0 1];
sequence = zeros(length(dataSoft),1);

for seq = 1 : length(dataSoft)
    tem = xor(scrambler(1),scrambler(4));
    scrambler(1:end-1) = scrambler(2:end);
    scrambler(end) = tem;
    sequence(seq) = tem;
end

descrambledData = xor(dataSoft,sequence);

% err_f = sum(abs((int8(descrambledData)-VHTData)))
% err_f_refe = sum(abs((int8(receive_descramble)-VHTData)))
% BER = err_f_refe/length(VHTData)

% DSP descrambling sequence creation
sequence = int8(sequence);
sequence_stuff_zero = [sequence;zeros(24,1)];
sequence_zero_reshape = reshape(sequence_stuff_zero,32,numel(sequence_stuff_zero)/32);
sequence_=sequence_zero_reshape';
sequence_deci=bi2de(sequence_);

% Get result
exclusive_or = xor(dataSoft_stuff_zero,sequence_stuff_zero);
sexclusive_orreshape = reshape(exclusive_or,32,numel(exclusive_or)/32);
sexclusive_orreshape_=uint8(sexclusive_orreshape');
reverse_result_wifi_=bi2de(sexclusive_orreshape_);

% isequal(reference_scramble,reverse_result_wifi_)
% descrambledData_tst = [descrambledData;zeros(24,1)];
% isequal(descrambledData_tst,exclusive_or)
