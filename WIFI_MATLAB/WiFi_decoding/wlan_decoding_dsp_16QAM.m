%% Initialization
% clear;
% clc;
% close all;

Nbpscs = 4; %% 8: 256QAM   4: 16QAM
R = 3/4;
% Coded bits per symbol (52 subcarrier * 8bits(256 QAM))
Ncbps = 52 * Nbpscs;
%% Data Convertion
% VarName1_reshape = reshape(VarName1,2,length(VarName1)/2);
% imag_data = VarName1_reshape(2,:);
% real_data = VarName1_reshape(1,:);
% waveform_data = real_data + 1j*imag_data;
%                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               
rxSig = synchFrame;  %% synchFrame  WIFI_Waveform RX_Data   wholeWaveform_scaling

%% Create waveform 
% tgacChan = wlanTGacChannel('SampleRate',20e6,'ChannelBandwidth','CBW20','DelayProfile','Model-C');
% 
% chNoise = comm.AWGNChannel('NoiseMethod','Signal to noise ratio (SNR)',...
%     'SNR',35);
% pfOffset = comm.PhaseFrequencyOffset('SampleRate',20e6,'FrequencyOffsetSource','Input port');
% % rxSig = chNoise(tgacChan(wholeWaveform));  %% wholeWaveform_scaling  wholeWaveform
% % rxSig = chNoise(wholeWaveform_scaling);
% rxSig = pfOffset(rxSig,1);
% rxSig = tgacChan(chNoise(pfOffset(rxSig, 200)));

%% DSP RF Data Preparation (NOT FOR THIS MATLAB SIMULATION)
% NOTE: The rxSig is from the synchronized frame sent by USRP
% Therefore, it is already scaled.
rf_data_int16 = int16(rxSig);
rf_data_int16_real = real(rf_data_int16);
rf_data_int16_imag = imag(rf_data_int16);
rf_data_ = [rf_data_int16_imag';rf_data_int16_real'];
rf_data_combine = reshape(rf_data_,length(rxSig)*2,1);

% Just put the following data into the header of DSP project to be the
% input data. (NOTE: it has four same frames)
rf_data_combine_4frames = repmat(rf_data_combine,[4,1]);
rf_data_combine_4frames_scale = rf_data_combine_4frames;
max___rx = max(abs(rf_data_combine_4frames_scale));

% Add offset to this combined four frame.
% Note:In order to make the synchronization result not to be zero
rf_data_combine_4frames_offset = [rf_data_combine_4frames(1560:end);rf_data_combine_4frames(1:1559)];


%% Remove CP(VHT-LTF)
vhtltf_part = rxSig(641:720,:);
vhtltf_cp_remove = vhtltf_part(17:end);

%% FFT( VHT-LTF )
% RF data
% IMPORTANT NOTE: The following way to calculate fft result is different from
% the standard way stated in the document of IEEE 802.11ac. The standard
% way should also be multiplied by sqrt(56)/64

% Without scaling sqrt(56)/64 (when doing channel compensation, because both denominator and
%numerator contain this term; therefore, it will be cancelled out)

% Another NOTE: for decoding, it should use fft, then fftshif. Otherwise,
% it will not conform with the standard one. (it would with 32 left-right shift)
fftc_dynamic_scaling = 2^3;
dsp_test_vhtltf = fft(vhtltf_cp_remove);
dsp_test_vhtltf_fftshift = double(fftshift(dsp_test_vhtltf,1));
dsp_test_vhtltf_fftshift = dsp_test_vhtltf_fftshift / fftc_dynamic_scaling;

% To check the dynamic scaling of FFTC in DSP
dsp_test_vhtltf_fftshift_1 =  int32(dsp_test_vhtltf_fftshift ) ;
dsp_test_vhtltf_fftshift_2 =  int32(dsp_test_vhtltf_fftshift * 2) ;
dsp_test_vhtltf_fftshift_1_2 =  int32(dsp_test_vhtltf_fftshift / 2) ;
dsp_test_vhtltf_fftshift_1_4 =  int32(dsp_test_vhtltf_fftshift / 4) ;
dsp_test_vhtltf_fftshift_1_8 =  int32(dsp_test_vhtltf_fftshift / 8) ;
 
%% Remove zero padding ( VHT-LTF )
vht_llf_af_pilots_scaling = [dsp_test_vhtltf_fftshift(5:32);dsp_test_vhtltf_fftshift(34:61)];

%% Channel Estimation ( VHT-LTF )
seqLong_VHT = [1;1;1;1;-1;-1;1;1;-1;1;-1;1;1;1;1;1;1;
    -1;-1;1;1;-1;1;-1;1;1;1;1;1;-1;-1;1;
    1;-1;1;-1;1;-1;-1;-1;-1;-1;1;1;-1;-1;1;-1;
    1;-1;1;1;1;1;-1;-1];
% Channel estimation with scaling factor.(DSP data)
H_VHT_scaling = vht_llf_af_pilots_scaling .*seqLong_VHT;

%% Remove CP(data)
data_part = rxSig(801:end,:);
res_data_bf_CP = reshape(data_part,[80,Nsym]);
data_cp_remove= res_data_bf_CP(17:80,:);


%% FFT (DATA FIELD)
% RF data
% IMPORTANT NOTE: The following way to calculate fft result is different from
% the standard way stated in the document of IEEE 802.11ac. The standard
% way should also be multiplied by sqrt(56)/64

% Without scaling sqrt(56)/64 (when doing channel compensation, because both denominator and
%numerator contain this term; therefore, it will be cancelled out)
dsp_test_data = fft(data_cp_remove);
dsp_test_fftshift = fftshift(dsp_test_data,1);
% dsp_test_fftshift = dsp_test_fftshift / fftc_dynamic_scaling;


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
zeropadding_remove = [dsp_test_fftshift(5:32,:);dsp_test_fftshift(34:61,:)];

%% Channel Compensation 
aft_estcom = zeros(56,Nsym);
% NOTE?The following way can not be used for DSP, because of the scaling of the
% result. Numerator and denominator have the same degree which results in
% the 0 degree result
for estnn = 1:Nsym
    aft_estcom(:,estnn) =  (zeropadding_remove(:,estnn))./(H_VHT_scaling);
end

% Test new way for dsp (Instead of dividing the denominator, multiply it
% with sqrt(1/170), which is the coefficient in 256QAM de-modulation.
compensated_real = zeros(56,Nsym);
compensated_imag = zeros(56,Nsym);
test_vht_real = real(H_VHT_scaling ) ;
test_vht_imag = imag(H_VHT_scaling ) ;
% The reason why 2^14 is chosen is out of experiment( appropriate scaling )
% divided by 64 is to divide the square of fftc dynamic scaling 8,
% 8*8=64

% NOTE: In normal way, nom should be divided by compensated_real and
% compensated_imag. Since division operation is too slow in DSP, and also
% because we are using SOFT de-mapping, the size of the value
% does not matter. However,in de-mapping algorithm, there is
% non-linear factor A, annd A has to be used in Addition operation. If without
% taking nom into account at all, the lack of nom would make A disproportional
% to the size of compensated_real and compensated_imag, which would lead to
% the failur of de-mapping. To avoid that, we just need to multiply nom
% with A to make them stay proportional to each other.

self_scaling = 2^13;

nom = (test_vht_real .* test_vht_real + test_vht_imag.* test_vht_imag)/ self_scaling;
nom_max = int16(max(nom));

for estnn = 1:Nsym
    aaa = real(zeropadding_remove(:,estnn)) .* test_vht_real;
    bbb = imag(zeropadding_remove(:,estnn)) .* test_vht_imag;
    ccc = imag(zeropadding_remove(:,estnn)) .* test_vht_real; % bc
    ddd = real(zeropadding_remove(:,estnn)) .* test_vht_imag; % ad
    compensated_real(:,estnn) = ((aaa + bbb)) / self_scaling;
    compensated_imag(:,estnn) = ((ccc - ddd)) / self_scaling;
end
compenstated_whole = compensated_real + i*compensated_imag;
compensated_real_max = int16(max(compensated_real));

compensated_real_reference = compensated_real/8;
compensated_imag_reference = compensated_imag / 8;

%% Residual Frequency Offset Estimation
% The following pilot is without nom, just numeritor. However, in fact, we
% just need numeriter (real part and imaginary part) to find the angle
pilot_noscale = [zeropadding_remove(8,:);zeropadding_remove(22,:);zeropadding_remove(35,:);zeropadding_remove(49,:)];
pilot_scale = [compenstated_whole(8,:);compenstated_whole(22,:);compenstated_whole(35,:);compenstated_whole(49,:)];

H_pilot = [H_VHT_scaling(8);H_VHT_scaling(22);H_VHT_scaling(35);H_VHT_scaling(49)];

residual_offset_dsp = zeros(Nsym,1);
residual_offset_test = zeros(Nsym,1);

% DSP applicable: The multiplication between conjugate of VHT and data
% pilot has already been calculated in channel compensation stage. It
% happens that it does not require to divide with denominator. All we need
% to do just add up the numerators which are the individual pilots
% multiplied with conjugated vht 
pilot_scale_polo = zeros(4,Nsym);
for resnn =1:Nsym
    pilo_ = pv(resnn+4)*pil(:,mod((resnn-1),4)+1);
    pilot_scale_polo(:,resnn) = pilot_scale(:,resnn).* pilo_;
end

pilot_scale_sum = sum(pilot_scale_polo,1);

% calculate angle
real_rcfo = real(pilot_scale_sum);
imag_rcfo = imag(pilot_scale_sum);

angle_RCFO_ = atan(imag_rcfo./real_rcfo);
angle_RCFO__ = angle_RCFO_'; % confirmed, it agrees with the result from the matlab code.


%% Residual Frequency Offset compensation
cos_rcfo = cos(angle_RCFO_);
sin_rcfo = sin(angle_RCFO_);

compensated_real_dsp = zeros(56,Nsym); 
compensated_imag_dsp = zeros(56,Nsym); 

for estnn = 1:Nsym
    compensated_real_dsp(:,estnn) = compensated_real(:,estnn).* cos_rcfo(estnn) + compensated_imag(:,estnn).* sin_rcfo(estnn);
    compensated_imag_dsp(:,estnn) = compensated_imag(:,estnn).* cos_rcfo(estnn) - compensated_real(:,estnn).* sin_rcfo(estnn);
end

% compensated_real_dsp_refer = compensated_real_dsp /8;
% compensated_imag_dsp_refer = compensated_imag_dsp / 8;
% pilot_remove_real_test_refe =[compensated_real_dsp_refer(1:7,:);compensated_real_dsp_refer(9:21,:);compensated_real_dsp_refer(23:34,:);compensated_real_dsp_refer(36:48,:);compensated_real_dsp_refer(50:56,:)];
% pilot_remove_imag_test_refe =[compensated_imag_dsp_refer(1:7,:);compensated_imag_dsp_refer(9:21,:);compensated_imag_dsp_refer(23:34,:);compensated_imag_dsp_refer(36:48,:);compensated_imag_dsp_refer(50:56,:)];

% compensated_imag_dsp = compensated_imag;
% compensated_real_dsp = compensated_real;

%% Remove pilots
% NOT FOR DSP (just for reference)
pilot_remove =[aft_estcom(1:7,:);aft_estcom(9:21,:);aft_estcom(23:34,:);aft_estcom(36:48,:);aft_estcom(50:56,:)];

% DSP test
pilot_remove_real_test =[compensated_real_dsp(1:7,:);compensated_real_dsp(9:21,:);compensated_real_dsp(23:34,:);compensated_real_dsp(36:48,:);compensated_real_dsp(50:56,:)];
pilot_remove_imag_test =[compensated_imag_dsp(1:7,:);compensated_imag_dsp(9:21,:);compensated_imag_dsp(23:34,:);compensated_imag_dsp(36:48,:);compensated_imag_dsp(50:56,:)];
pilot_remove_nom = [nom(1:7,:);nom(9:21,:);nom(23:34,:);nom(36:48,:);nom(50:56,:)];
pilot_remove_nom_1_2 = pilot_remove_nom;
% The following A_scaling will be used in De-mapping to make the divion
% converted to multiplication possible without affecting the proportion
% between A and compensated_real, compensated_imag
% NOTE: If Encoding project does not multiple A when doinng mapping, then
% decoding process without A=sqrt(1/170) (treat A as 1) can still get correct result
A_scaling  =  sqrt(1/170) * pilot_remove_nom_1_2;
% A_scaling  =  1 * pilot_remove_nom_1_2;

% This result should be the same as pilot_remove (just for reference)
pilot_remove_test = (pilot_remove_real_test+1i*pilot_remove_imag_test)./pilot_remove_nom ;
pilot_remove_test_one = reshape(pilot_remove_test,numel(pilot_remove_test),1);


% % test (delete later)
% RX_Data_r = reshape(RX_Data,52,53);
% RX_Data_r_de = RX_Data_r./pilot_remove_nom;
% RX_Data_r_de_one = reshape(RX_Data_r_de,numel(RX_Data_r_de),1);
sss = 8;

% plot the result of compensated data
figure(sss)
xlim([-1.5 1.5])
ylim([-1.5 1.5])
% pilot_remove_test_one_ = pilot_remove_test_one;
sz = 25;
% scatter(real(pilot_remove_test(:,sss)),imag(pilot_remove_test(:,sss)),sz,'*') %pilot_remove_test_one  RX_Data_r_de_one
scatter(real(pilot_remove_test_one),imag(pilot_remove_test_one),sz,'*') %pilot_remove_test_one  RX_Data_r_de_one
xlabel('In-phase')
ylabel('Quadrature')
title('Received 256 QAM Constellation')

%% Constellation Demapper(algorithm)
llrAggregate = zeros(Ncbps, Nsym);
llr = zeros(8,1);

qamBits= zeros(Nbpscs,52);
for oo = 1  :  Nsym
    for ll = 1: 52
        A = A_scaling(ll);
        rr = pilot_remove_real_test(ll,oo);
        im = pilot_remove_imag_test(ll,oo);
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
after_puncturing = zeros(8424*2,1);
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

% NOTE: The following equation checks whether the size exceeds the range of
% VCP2. Very important step to know how to set the scaling of BM.
BM_scale = BM/2^8; 
[maxi,ind] = max(abs(BM_scale));

BM_scale_int8 = int8(BM_scale);

BM_scale_reshape = reshape(BM_scale_int8,624,Nsym);

%% BCC Decoder
% WIFI decoder function (for reference to dataSoft)
decodedData = wlanBCCDecode(puncturing_reshape,1/2,'soft');

% Viterbi decoder setting
trellis = poly2trellis(7,[133 171]);
tbl = 32;
rate = 1/2;
% Viterbi decode the demodulated data
% dataHard = vitdec(rxDataHard,trellis,tbl,'cont','hard');
dataSoft = vitdec(puncturing_reshape,trellis,tbl,'trunc','unquant');

% Below is for DSP result confirmation (Decisions of VCP2)
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
err_f = sum(abs((double(descrambledData)-VHTData)))

% DSP descrambling sequence creation ( For DSP: put into header)
% NOTE: Add extra zero to make it aligned with 32bit
sequence = int8(sequence);
sequence_stuff_zero = [sequence;zeros(24,1)];
sequence_zero_reshape = reshape(sequence_stuff_zero,32,numel(sequence_stuff_zero)/32);
sequence_=sequence_zero_reshape';
sequence_deci=bi2de(sequence_);

% To verify the final result of DSP (FOR DSP)
exclusive_or = xor(dataSoft_stuff_zero,sequence_stuff_zero);
sexclusive_orreshape = reshape(exclusive_or,32,numel(exclusive_or)/32);
sexclusive_orreshape_=uint8(sexclusive_orreshape');
reverse_result_wifi_=bi2de(sexclusive_orreshape_);
