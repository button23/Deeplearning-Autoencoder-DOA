function pilot_remove = IPS_dataDecoding(rxSig, cfgVHT, H_VHT)
%% Remove CP(data)
ind_Data = wlanFieldIndices(cfgVHT,'VHT-Data');

if cfgVHT.MCS ==8
    Nsym =27;
else
    Nsym = 53;
end

data_part = rxSig(ind_Data(1):ind_Data(2), :);
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

%% Remove zero padding
zeropadding_remove_origin = [vhtdata_fft(5:32,:);vhtdata_fft(34:61,:)];

%% Channel Compensation and Residual Frequency Offset compensation
aft_estcom = zeros(56,Nsym);

for estnn = 1:Nsym
    %    aft_estcom(:,estnn) =  zeropadding_remove_origin(:,estnn)./(H_VHT) * exp(-1i*residual_offset(estnn));
    aft_estcom(:,estnn) =  (zeropadding_remove_origin(:,estnn))./(H_VHT) ;   %%  This one can be used to observe the CFO
end

%% Remove pilots
pilot_remove =[aft_estcom(1:7,:);aft_estcom(9:21,:);aft_estcom(23:34,:);aft_estcom(36:48,:);aft_estcom(50:56,:)];

% plot the result of compensated data
figure(6)
xlim([-1.5 1.5])
ylim([-1.5 1.5])
pilot_remove_ =(reshape(pilot_remove,[52*Nsym,1]));
sz = 25;
scatter(real(pilot_remove_),imag(pilot_remove_),sz,'*')
xlabel('In-phase')
ylabel('Quadrature')
title('Received 16 QAM Constellation')

end