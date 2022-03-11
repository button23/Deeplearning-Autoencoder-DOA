%% Initialization
% clear;
% clc;
close all;

%% Data Convertion
rxSig = synchFrame;  %% synchFrame frame_sync   rxSig WIFI_Waveform  RX_Data

%% wlan waveform
cfgVHT = wlanVHTConfig('ChannelBandwidth','CBW20','MCS',4);
% cfgVHT
% Coded bits per single carrier for each spatial steam
Nbpscs = 4; %% 8: 256QAM   4: 16QAM
R = 3/4;
% Coded bits per symbol (52 subcarrier * 8bits(256 QAM))
Ncbps = 52 * Nbpscs;
% Data bits per OFDM symbol (416 * 3/4 = 312 bits)
Ndbps = Ncbps * R;
% Used for scaling Matlab data in order to move to DSP
scalingFactor = 2^13;
% Default value is 1024
cfgVHT.APEPLength = 1024;
% Create random bit stream
rng(25)
PSDU = randi([0 1],cfgVHT.PSDULength*8,1);
% PSDU = [1;0;0;1];
% PSDU = repmat(PSDU,2100,1);
% Convert PSDU into ppdu
wholeWaveform = wlanWaveformGenerator(PSDU,cfgVHT);
wholeWaveform_scaling = 2^13 * wholeWaveform;
% wholeWaveform_scaling(641:720) = wholeWaveform_scaling(641:720)/(64/sqrt(56));
% wholeWaveform_scaling(801:end) = wholeWaveform_scaling(801:end)/(64/sqrt(56));

% [startOffset,meric] = wlanSymbolTimingEstimate(wholeWaveform,'CBW20'); %% rxSig  RX_Data WIFI_Waveform
% % RX_Data = RX_Data(startOffset+1:startOffset+2960);
% figure(22)
% plot(meric)

ind_Data = wlanFieldIndices(cfgVHT,'VHT-Data');
ind_VHTSIGA = wlanFieldIndices(cfgVHT,'VHT-SIG-A');
ind_VHTSIGB = wlanFieldIndices(cfgVHT,'VHT-SIG-B');

ind_VHTLTF = wlanFieldIndices(cfgVHT,'VHT-LTF');
ind_VHTSTF = wlanFieldIndices(cfgVHT,'VHT-STF');

ind_STF = wlanFieldIndices(cfgVHT,'L-STF');
ind_LTF = wlanFieldIndices(cfgVHT,'L-LTF');
ind_LSIG = wlanFieldIndices(cfgVHT,'L-SIG');

data_part_transmit = wholeWaveform(ind_Data(1):ind_Data(2),:);

% PSDU = int8(PSDU);
sigDATA = wlanVHTData(PSDU,cfgVHT);  % 93 default     2160= 801:2960

%% Preamble
sigLSTF = wlanLSTF(cfgVHT);
sigLLTF = wlanLLTF(cfgVHT);
[sigLSIG, bits_sigLSIG] = wlanLSIG(cfgVHT);
[sigA,bits_sigA] = wlanVHTSIGA(cfgVHT);
[sigB,bits_sigB] = wlanVHTSIGB(cfgVHT);
sigVHTSTF = wlanVHTSTF(cfgVHT);
sigVHTLTF = wlanVHTLTF(cfgVHT);

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
% Insert padding bits and service field into data rfield
VHTData = [VHTData_Service;PSDU;zeros(Npad,1);zeros(6,1)];

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
% rxSig = pfOffset(rxSig,250);
% rxSig = tgacChan(chNoise(pfOffset(wholeWaveform, 28000)));
% rxSig = pfOffset(wholeWaveform, 300);

% rxSig = WIFI_Waveform;  %% WIFI_Waveform  wholeWaveform

%% Data recover (WLAN function)  ind_VHTLTF
% ltfDemodSig = wlanVHTLTFDemodulate(synchFrame(ind_VHTLTF(1):ind_VHTLTF(2),:), cfgVHT);
% chEst = wlanVHTLTFChannelEstimate(ltfDemodSig,cfgVHT);
% recBits = wlanVHTDataRecover(synchFrame(801:end),chEst,0,cfgVHT);
% numErrs = biterr(PSDU,recBits)

% figure(40)
% plot(abs(chEst))
% hold on
% xlabel('Subcarrier')
% ylabel('Magnitude')

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
rxSig(ind_LTF(1)+32:ind_LTF(2)) = rxSig(ind_LTF(1)+32:ind_LTF(2)) .* exp(-1i * mmC' * phase_STF);

%% Fine CFO Estimation
Ts = 1/(2*10^7);
ltf_part = rxSig(ind_LTF(1):ind_LTF(2),:);
first_LTF = ltf_part(33:96);
second_LTF = ltf_part(97:160);

phase_LTF = 1/64 * angle(first_LTF'*second_LTF);
CFO_LTF = phase_LTF/(2 * pi * Ts)
% % Total carrier frequency offset phase combined with coarse CFO estimation and
% % Fine CFO estimation
total_phase = phase_LTF + phase_STF;
total_cfo = CFO_LTF + CFO_STF

%% Wlan Funciton to estimate CFO
% WIFI function
% foffset1 = wlanCoarseCFOEstimate(stf_part,'CBW20')
% rxSig= pfOffset(rxSig,-foffset1);
% ltf_part = rxSig(ind_LTF(1):ind_LTF(2),:);
% foffset2 = wlanFineCFOEstimate(ltf_part,'CBW20')
% rxSig= pfOffset(rxSig,-foffset2);

%% Fine CFO Correlation
mmF = 128:128+ind_Data(2)-ind_LTF(2)-1;
% mmF = 128:ind_Data(2);

mmF = double(mmF');

% If do Fine Compensation at the beginning, the error rate is less and the
% constallation map is  neater.
rxSig(ind_LSIG(1):end) = rxSig(ind_LSIG(1):end) .* exp(-j * mmF * total_phase);
% rxSig(ind_LTF(2)+1:end) = rxSig(ind_LTF(2)+1:end) .* exp(-j * mmF * phase_STF);
% rxSig(ind_LTF(2)+1:end) = rxSig(ind_LTF(2)+1:end) .* exp(-j * mmF * phase_LTF);

%% L-LTF FFT
rxSig_lega_ltf = rxSig(ind_LTF(1):ind_LTF(2));
rxSig_lega_ltf_ = reshape(rxSig_lega_ltf,80,2);
rxSig_lega_ltf_rcp = rxSig_lega_ltf_(17:end,:);
rxSig_lega_ltf_fft = fft(rxSig_lega_ltf_rcp);
rxSig_lega_ltf_fftshift = double(fftshift(rxSig_lega_ltf_fft,1))*1/64 *sqrt(52);
% Remove zero padding ( L-LTF )
l_llf_fft = [rxSig_lega_ltf_fftshift(7:32,:);rxSig_lega_ltf_fftshift(34:59,:)];

%% Channel Estimation (L-LTF)
seqLong_L = [1;1;-1;-1;1;1;-1;1;-1;1;1;1;1;1;1;
    -1;-1;1;1;-1;1;-1;1;1;1;1;1;-1;-1;1;
    1;-1;1;-1;1;-1;-1;-1;-1;-1;1;1;-1;-1;1;-1;
    1;-1;1;1;1;1];
H_L = (l_llf_fft(:,1)+l_llf_fft(:,2))/2 .* seqLong_L;

%% L-SIG FFT
rxSig_lega_sig = rxSig(ind_LSIG(1):ind_LSIG(2));
rxSig_lega_sig_rcp = rxSig_lega_sig(17:end);
rxSig_lega_sig_fft = fft(rxSig_lega_sig_rcp);
rxSig_lega_sig_fftshift = double(fftshift(rxSig_lega_sig_fft,1))*1/64 *sqrt(52);
% remove zero padding
lega_sig = [rxSig_lega_sig_fftshift(7:32);rxSig_lega_sig_fftshift(34:59)];

%% VHT-SIG-A FFT
rxSig_vht_sig = rxSig(ind_VHTSIGA(1):ind_VHTSIGA(2));
rxSig_vht_sig_ = reshape(rxSig_vht_sig,80,2);
rxSig_vht_sig_rcp = rxSig_vht_sig_(17:end,:);
rxSig_vht_sig_fft = fft(rxSig_vht_sig_rcp);
rxSig_vht_sig_fftshift = double(fftshift(rxSig_vht_sig_fft,1))*1/64 *sqrt(52);
% remove zero padding
lega_sig_a = [rxSig_vht_sig_fftshift(7:32,:);rxSig_vht_sig_fftshift(34:59,:)];

%% VHT-STF FFT
rxSig_vht_stf = rxSig(ind_VHTSTF(1):ind_VHTSTF(2));
rxSig_vht_stf_rcp = rxSig_vht_stf(17:end,:);
rxSig_vht_stf_fft = fft(rxSig_vht_stf_rcp);
rxSig_vht_stf_fftshift = double(fftshift(rxSig_vht_stf_fft,1))*1/64 *sqrt(12);
% remove zero padding
vht_stf = [rxSig_vht_stf_fftshift(5:32);rxSig_vht_stf_fftshift(34:61)];

%% VHT-SIG-B FFT
rxSig_vht_sigb = rxSig(ind_VHTSIGB(1):ind_VHTSIGB(2));
rxSig_vht_sigb_rcp = rxSig_vht_sigb(17:end,:);
rxSig_vht_sigb_fft = fft(rxSig_vht_sigb_rcp);
rxSig_vht_sigb_fftshift = double(fftshift(rxSig_vht_sigb_fft,1))*1/64 *sqrt(56);
% remove zero padding
vht_sig_b = [rxSig_vht_sigb_fftshift(5:32);rxSig_vht_sigb_fftshift(34:61)];

%% VHT-LTF FFT
vhtltf_part = rxSig(ind_VHTLTF(1):ind_VHTLTF(2),:);
cp_remove_vhtltf = vhtltf_part(17:end);
% k2 = 1:1:64;
% k2 = k2';
% vht_llf_fft = zeros(64,1);
%
% for kk = 1:64
%     yy =  sqrt(56) *cp_remove_vhtltf(:) .* exp(-j*2*pi*(kk-33)*(k2-1)/64);
%     vht_llf_fft(kk) = 1/64 *  sum(yy) ;
% end

% just for reference
dsp_test_vhtltf = fft(cp_remove_vhtltf);
dsp_test_vhtltf_fftshift = double(fftshift(dsp_test_vhtltf,1))*1/64 *sqrt(56);
% Remove zero padding ( VHT-LTF )
vht_llf_af_pilots = [dsp_test_vhtltf_fftshift(5:32);dsp_test_vhtltf_fftshift(34:61)];

%% Channel Estimation ( VHT-LTF )
seqLong_VHT = [1;1;1;1;-1;-1;1;1;-1;1;-1;1;1;1;1;1;1;
    -1;-1;1;1;-1;1;-1;1;1;1;1;1;-1;-1;1;
    1;-1;1;-1;1;-1;-1;-1;-1;-1;1;1;-1;-1;1;-1;
    1;-1;1;1;1;1;-1;-1];
H_VHT = vht_llf_af_pilots .* seqLong_VHT;

%% Remove CP(data)
data_part = rxSig(ind_Data(1):ind_Data(2),:);
res_data_bf_CP = reshape(data_part,[80,Nsym]);
data_cp_remove= res_data_bf_CP(17:80,:);

%% FFT (data)
k2 = 1:1:64;
k2 = k2';
vhtdata_fft = zeros(64,Nsym);
for fftnn = 1:Nsym
    for kk = 1:64
        yy =  sqrt(56) *data_cp_remove(:,fftnn) .* exp(-j*2*pi*(kk-33)*(k2-1)/64);
        vhtdata_fft(kk,fftnn) = 1/64 *  sum(yy) ;
    end
end

% just for reference
dsp_test_data = fft(data_cp_remove);
dsp_test_fftshift = fftshift(dsp_test_data,1) * sqrt(56) * 1/64;


%% Remove zero padding
zeropadding_remove_origin = [vhtdata_fft(5:32,:);vhtdata_fft(34:61,:)];
residual_offset = zeros(Nsym,1);
pilot_data_origin = zeros(4,1);
residualPhase = 0;
Tu = 4 * 10^(-6); % symbol time
fc = 2.4 * 10^9;
fs = 20 * 10^6;
totalPhase = total_phase/(2*pi*fc*Tu);
weightValueAccu = 0;

% Pilots in VHT-LTF     
H_pilot = [H_VHT(8);H_VHT(22);H_VHT(35);H_VHT(49)];
%     H_L_pilot = [H_L(6);H_L(20);H_L(33);H_L(47)];
%     Pilot_four = [1,1,1,-1];
% lega_sig_pilot = [lega_sig_cor(6);lega_sig_cor(20);lega_sig_cor(33);lega_sig_cor(47)];
%    rcfSig = angle(lega_sig_pilot'*H_L_pilot);
offIn1 = -28:-1;
offIn2 = 1:28;
offIn = [offIn1,offIn2];

pilot_data_middle = H_pilot;
for sfoc =1:Nsym
    %% Step 1: Correct Sampling Frequency Offset after FFT
    %The estimated phase alpha is also used to find an intial estimate epxilong of the frequency
    totalPhase =totalPhase + residualPhase;
%         totalPhase =totalPhase;

%     phasss = total_phase + totalPhase;
%     % lega_sig_cor =  lega_sig.*exp(2 * pi * 80 / 64 * initial_e * offIn');
%     zeropadding_remove_origin(:,sfoc) =  zeropadding_remove_origin(:,sfoc).*exp(2 * pi * (7 + sfoc) * 80 / 64 * totalPhase * offIn');
    
    %% Step 2: Residual Frequency Offset Estimation
    pilot_data_origin = [zeropadding_remove_origin(8,sfoc);zeropadding_remove_origin(22,sfoc);zeropadding_remove_origin(35,sfoc);zeropadding_remove_origin(49,sfoc)];
    % Multiplication by its corresponding pseudo-random value to remove the
    % randomness
    pilo_ = pv(sfoc+4)*pil(:,mod((sfoc-1),4)+1);
    % Save the RFO for all the OFDM symbols, used in channel compensation
    % step
    residual_offset(sfoc) = angle((H_pilot.*pilo_)'*pilot_data_origin);
    
    %% Step 3: Fine Estimation of Frequency Offset
    % No filtering
    weight_updata = [pilot_data_middle,pilot_data_origin];
    weightValue = weight_updata(:,1)'*weight_updata(:,2);
    % Save pilot for next use
    pilot_data_middle = pilot_data_origin;
    weightValueAccu = weightValueAccu + weightValue;
    residualPhase= angle(weightValueAccu) / (2 * pi * Tu * fc);
     
    % Filtering
    % for nnns = 1:length(weightAdd)
    %     weightAdd(nnns) = sum(weight((4*(nnns-1)+1):4*nnns));
    % end

    % pp = 1/32;
    % weightFilter = pp
end

%% Remove pilots
pilot_remove =[zeropadding_remove_origin(1:7,:);zeropadding_remove_origin(9:21,:);zeropadding_remove_origin(23:34,:);zeropadding_remove_origin(36:48,:);zeropadding_remove_origin(50:56,:)];
H_VHT_rem_pilot = [H_VHT(1:7);H_VHT(9:21);H_VHT(23:34);H_VHT(36:48);H_VHT(50:56)];

%% Channel Compensation and Residual Frequency Offset compensation
aft_estcom = zeros(length(H_VHT_rem_pilot),Nsym);

mat_cos = cos(residual_offset);
mat_sin = sin(residual_offset);

for estnn = 1:Nsym
    aft_estcom(:,estnn) =  pilot_remove(:,estnn)./(H_VHT_rem_pilot) * exp(-1i*residual_offset(estnn));
%         aft_estcom(:,estnn) =  (pilot_remove(:,estnn))./(H_VHT_rem_pilot) ;   %%  This one can be used to observe the CFO
end


%
% aft_estcom_real = real(aft_estcom);
% aft_estcom_imag = imag(aft_estcom);
% aft_estcom_real_ = zeros(56,Nsym);
% aft_estcom_imag_ = zeros(56,Nsym);
%
% % % %  Residual offset compnesation
% for estnn = 1:Nsym
%     aft_estcom_real_(:,estnn) =  aft_estcom_real(:,estnn) * mat_cos(estnn) - aft_estcom_imag(:, estnn) * mat_sin(estnn);
%     aft_estcom_imag_(:,estnn) =  aft_estcom_real(:,estnn) * mat_sin(estnn) + aft_estcom_imag(:, estnn) * mat_cos(estnn);
% end
%
% aft_estcom = aft_estcom_real_ + i*aft_estcom_imag_;


%% plot the result of compensated data
figure(6)
lslsls = 20;

xlim([-1.5 1.5])
ylim([-1.5 1.5])
pilot_remove_ =(reshape(aft_estcom,[52*Nsym,1]));
sz = 25;
% scatter(real(pilot_remove(:,lslsls)),imag(pilot_remove(:,lslsls)),sz,'*')
scatter(real(pilot_remove_),imag(pilot_remove_),sz,'*')
xlabel('In-phase')
ylabel('Quadrature')
title('Received 16 QAM Constellation')

%% Constellation Demapper(algorithm)
llrAggregate = zeros(Ncbps, Nsym);
llr = zeros(8,1);
A = sqrt(1/170);

qamBits= zeros(Nbpscs,52);
for oo = 1  :  Nsym
    for ll = 1: 52
        rr = real(pilot_remove(ll,oo));
        im = imag(pilot_remove(ll,oo));
        llr(1) = rr ;
        llr(2) = - abs(llr(1)) + 8 * A;
        llr(3) = - abs( llr(2) ) + 4 * A;
        llr(4) = - abs( llr(3)) + 2 * A;
        
        llr(5) = im ;
        llr(6) = - abs(llr(5)) + 8 * A;
        llr(7) = - abs( llr(6)) + 4 * A;
        llr(8) = - abs( llr(7)) + 2 * A;
        
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

%% BCC Decoder
% WIFI decoder function
decodedData = wlanBCCDecode(puncturing_reshape,1/2,'soft');

% Viterbi decoder setting
trellis = poly2trellis(7,[133 171]);
tbl = 32;
rate = 1/2;
% Viterbi decode the demodulated data
% dataHard = vitdec(rxDataHard,trellis,tbl,'cont','hard');
dataSoft = vitdec(puncturing_reshape,trellis,tbl,'trunc','unquant');

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

err_f = sum(abs((double(descrambledData)-VHTData)))
% err_f_refe = sum(abs((double(receive_descramble)-VHTData)))
% BER = err_f_refe/length(VHTData)
