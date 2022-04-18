% input_wlan = VarName17;  %% matlabdata  dspdata  VarName2 rf_data_combine_4frames_scale
% input_reshape = reshape(input_wlan,2,numel(input_wlan)/2);
% real_input = input_reshape(2,:).';
% imag_input = input_reshape(1,:).';
% RX_Data = real_input+i*imag_input;
% RX_Data_int32 = int32(RX_Data);

frameLength = 350;
framewithCfo = zeros(5040,frameLength);
pilot_remove_acc = zeros(2756,frameLength);
pilot_remove_cfo_acc = zeros(2756,frameLength);

residual_commu  = zeros(53,frameLength);

for recomu = 1:frameLength
    
    synchFrame = RX_Data((synchResult_offset + (recomu-1) * 5040):(synchResult_offset+5039 + (recomu-1) * 5040)); %% the same as  frame_sync which has different noise in it
    
    rxSig = synchFrame;  %%frame_sync  synchFrame rxSig WIFI_Waveform  RX_Data
    
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
    CFO_STF = phase_STF/(2*pi*Ts);
    
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
    CFO_LTF = phase_LTF/(2 * pi * Ts);
    % % Total carrier frequency offset phase combined with coarse CFO estimation and
    % % Fine CFO estimation
    total_phase = phase_LTF + phase_STF;
    total_cfo = CFO_LTF + CFO_STF;
    
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
    
    %% Correct Sampling Frequency Offset ?after FFT?
    offIn1 = -28:-1;
    offIn2 = 1:28;
    offIn = [offIn1,offIn2];
    
    %The estimated phase alpha is also used to find an intial estimate epxilong of the frequency
    initial_e = total_phase;
    % lega_sig_cor =  lega_sig.*exp(2 * pi * 80 / 64 * initial_e * offIn');
    % for symmm = 7:7+Nsym
    %     zeropadding_remove_origin =  zeropadding_remove_origin.*exp(2 * pi * symmm * 80 / 64 * initial_e * offIn');
    % end
    %% Correct Residual Carrier Frequency Offset
    H_pilot = [H_VHT(8);H_VHT(22);H_VHT(35);H_VHT(49)];
    H_L_pilot = [H_L(6);H_L(20);H_L(33);H_L(47)];
    
    Pilot_four = [1,1,1,-1];
    % lega_sig_pilot = [lega_sig_cor(6);lega_sig_cor(20);lega_sig_cor(33);lega_sig_cor(47)];
    % rcfSig = angle(lega_sig_pilot'*H_L_pilot);
    
    % lega_sig_cor
    
    %% Residual Frequency Offset Estimation
    pilot_data_origin = [zeropadding_remove_origin(8,:);zeropadding_remove_origin(22,:);zeropadding_remove_origin(35,:);zeropadding_remove_origin(49,:)];
    
    residual_offset = zeros(Nsym,1);
    residual_offset_1 = zeros(Nsym,1);
    
    
    for resnn =1:Nsym
        pilo_ = pv(resnn+4)*pil(:,mod((resnn-1),4)+1);
        residual_offset(resnn) = angle((H_pilot.*pilo_)'*pilot_data_origin(:,resnn));
    end
    
    %% Fine Estimation of Frequency Offset
    weight = zeros(53,1);
    weightAdd = zeros(int16(Nsym/4),1);
    pilot_data_origin_2 = [H_pilot,pilot_data_origin];
    for kkres = 1:Nsym
        weight(kkres) = pilot_data_origin_2(:,kkres)'*pilot_data_origin_2(:,kkres+1);
    end
    % for kkres = 1:Nsym
    %     weight(kkres) = pilot_data_origin_2(:,kkres+1)'*pilot_data_origin_2(:,kkres);
    % end
    
    %% No filtering
    offIn1 = -28:-1;
    offIn2 = 1:28;
    offIn = [offIn1,offIn2];
    
    Tu = 4 * 10^-6; % symbol time
    fc = 2.4 * 10^9;
    epusil = angle(weight) / (2 * pi * Tu * fc);
    
    for kknn = 1 : Nsym
        %     zeropadding_remove_origin(:,kknn) =  zeropadding_remove_origin(:,kknn).*exp(2 * pi * (kknn+7) * 80 / 64 * epusil(kknn) * offIn');
    end
    
    %% Filtering
    % for nnns = 1:length(weightAdd)
    %     weightAdd(nnns) = sum(weight((4*(nnns-1)+1):4*nnns));
    % end
    
    % pp = 1/32;
    % weightFilter = pp
    
    %% Channel Compensation and Residual Frequency Offset compensation
    aft_estcom = zeros(56,Nsym);
    aft_estcom_cfo = zeros(56,Nsym);
    
    
    mat_cos = cos(residual_offset);
    mat_sin = sin(residual_offset);
    
    for estnn = 1:Nsym
        aft_estcom_cfo(:,estnn) =  zeropadding_remove_origin(:,estnn)./(H_VHT) * exp(-1i*residual_offset(estnn));
        
        aft_estcom(:,estnn) =  (zeropadding_remove_origin(:,estnn))./(H_VHT) ;   %%  This one can be used to observe the CFO
    end
    
    
    residual_commu(recomu,1) = total_cfo;
    
    %% Remove pilots
    pilot_remove =[aft_estcom(1:7,:);aft_estcom(9:21,:);aft_estcom(23:34,:);aft_estcom(36:48,:);aft_estcom(50:56,:)];
    pilot_remove_ =(reshape(pilot_remove,[52*Nsym,1]));
    pilot_remove_acc(:,recomu) = pilot_remove_;
    
    pilot_remove_cfo =[aft_estcom_cfo(1:7,:);aft_estcom_cfo(9:21,:);aft_estcom_cfo(23:34,:);aft_estcom_cfo(36:48,:);aft_estcom_cfo(50:56,:)];
    pilot_remove_cfo_ =(reshape(pilot_remove_cfo,[52*Nsym,1]));
    pilot_remove_cfo_acc(:,recomu) = pilot_remove_cfo_;
end
ratttt = zeros(2756,349);
for ssns=1:349
    ratttt(:,ssns) = pilot_remove_cfo_acc(:,ssns+1)./pilot_remove_cfo_acc(:,ssns);
end
