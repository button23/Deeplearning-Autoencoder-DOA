%% Initialize
% close all
% clear
% clc

%% 수신 데이터 complex value 구성
for i=1:153600 %307200 76800 153600
    imag_enb(i,:) = VarName1(2*i-1,:);
    real_enb(i,:) = VarName1(2*i,:);
end

waveform = real_enb + 1j*imag_enb;

waveform = waveform .*2^-3/2868;

%% processing
% reall = real(waveform)./abs(real(waveform));
% imagg = imag(waveform)./abs(imag(waveform));
% summ = reall + imagg;
% indd = summ==0;
% test = waveform;
% test(indd) = -1 * test(indd);


%% Cell-Wide Settings
rmccfg.RC = 'R.7';
ncodewords = 1;
enb = lteRMCDL(rmccfg,ncodewords);

%% CFO Estimation and Compensation
cfo_offset = lteFrequencyOffset(enb,waveform); %waveform waveform
disp(['Frequency offset: ' num2str(cfo_offset) ' Hz'])

waveform1_nocfo = lteFrequencyCorrect(enb,waveform,cfo_offset);
cfo_o = lteFrequencyOffset(enb,waveform1_nocfo); %waveform waveform
disp(['Frequency offset: ' num2str(cfo_o) ' Hz'])
%% Synchronization
close all
[offset,corr] = lteDLFrameOffset(enb,waveform1_nocfo); %%txWaveform rxWaveform waveform

figure();plot(corr)
%/////////////////////////////////////////////////////////////////////////
fft_size = 1024; cp = 72; one_sym = fft_size + cp;
% offset = offset + 9; % CP is always 72, add 8 extra zeros to make it 80 for the first subframe.
offset = offset + (-10); % CP is always 72, add 8 extra zeros to make it 80 for the first subframe.
subframe_data = waveform1_nocfo(offset+1 :offset + one_sym *14);


rxGrid_fft  = [fft(subframe_data(cp+1:one_sym)), fft(subframe_data(one_sym+cp+1:one_sym*2)), ...
            fft(subframe_data(one_sym*2+cp+1:one_sym*3)), fft(subframe_data(one_sym*3+cp+1:one_sym*4)), ...
            fft(subframe_data(one_sym*4+cp+1:one_sym*5)), fft(subframe_data(one_sym*5+cp+1:one_sym*6)), ...
            fft(subframe_data(one_sym*6+cp+1:one_sym*7)), fft(subframe_data(one_sym*7+cp+1:one_sym*8)), ...
            fft(subframe_data(one_sym*8+cp+1:one_sym*9)), fft(subframe_data(one_sym*9+cp+1:one_sym*10)), ...
            fft(subframe_data(one_sym*10+cp+1:one_sym*11)), fft(subframe_data(one_sym*11+cp+1:one_sym*12)), ...
            fft(subframe_data(one_sym*12+cp+1:one_sym*13)), fft(subframe_data(one_sym*13+cp+1:one_sym*14))];
        
sym1 = [rxGrid_fft(end-299:end,1);rxGrid_fft(2:301,1)];
sym2 = [rxGrid_fft(end-299:end,2);rxGrid_fft(2:301,2)];
sym3 = [rxGrid_fft(end-299:end,3);rxGrid_fft(2:301,3)];
sym4 = [rxGrid_fft(end-299:end,4);rxGrid_fft(2:301,4)];
sym5 = [rxGrid_fft(end-299:end,5);rxGrid_fft(2:301,5)];
sym6 = [rxGrid_fft(end-299:end,6);rxGrid_fft(2:301,6)];
sym7 = [rxGrid_fft(end-299:end,7);rxGrid_fft(2:301,7)];
sym8 = [rxGrid_fft(end-299:end,8);rxGrid_fft(2:301,8)];
sym9 = [rxGrid_fft(end-299:end,9);rxGrid_fft(2:301,9)];
sym10 = [rxGrid_fft(end-299:end,10);rxGrid_fft(2:301,10)];
sym11 = [rxGrid_fft(end-299:end,11);rxGrid_fft(2:301,11)];
sym12 = [rxGrid_fft(end-299:end,12);rxGrid_fft(2:301,12)];
sym13 = [rxGrid_fft(end-299:end,13);rxGrid_fft(2:301,13)];
sym14 = [rxGrid_fft(end-299:end,14);rxGrid_fft(2:301,14)];

rxGrid = [sym1, sym2, sym3, sym4, sym5, sym6, sym7,...
        sym8, sym9, sym10, sym11, sym12, sym13, sym14];
        % subframe_data = [subframe_data(1023-7:1023); ...
%     subframe_data(1: one_sym *7); ...
%     subframe_data(one_sym *7 +1023-7:one_sym *7 +1023); ...
%     subframe_data(one_sym *7 + 1:end)];
% % OFDM Demodulation
% rxGrid = lteOFDMDemodulate(enb,subframe_data); %% rxWaveform_cfo_corrected subframe_data
% % for i =1:14
% %     % plot test
% %     scatterplot(rxGrid(:,i));
% % end
% % plot cell reference signal
scatterplot(rxGrid(1:6:end,1)); % symbol 1
scatterplot(rxGrid(4:6:end,5));% symbol 5
scatterplot(rxGrid(1:6:end,8));% symbol 8
scatterplot(rxGrid(4:6:end,12));% symbol 12

%% Channel Estimation(Coding)
% Extract CRS from real data
% enb.NSubframe = 0;
cellRsInd = lteCellRSIndices(enb);
crs_subfram = rxGrid(cellRsInd); % get crs from received signal

cellRsSym = lteCellRS(enb,0);

scatterplot(crs_subfram);

% h_est_coding = crs_subfram ./ cellRsSym; % channel H
% figure(29)
% % plot(abs(h_est_coding))
% scatter(real(h_est_coding),imag(h_est_coding))
% 
% 
% scatterplot(h_est_coding(1:5))
% % scatterplot(h_est_coding(101:200))
% scatterplot(h_est_coding(201:300))
% scatterplot(h_est_coding(301:400))


%% processing
reall = real(crs_subfram)./abs(real(crs_subfram));
imagg = imag(crs_subfram)./abs(imag(crs_subfram));
summ = reall + imagg;
indd = summ==0;
test = crs_subfram;
test(indd) = -1 * test(indd);

h_est_coding = test ./ cellRsSym; % channel H
% figure(29)
% plot(abs(h_est_coding))
% scatter(real(h_est_coding),imag(h_est_coding))
%%
% % Time averaging
% crs_subfram_res = reshape(crs_subfram,200,4);
% ave1 = (crs_subfram_res(:,1) + crs_subfram_res(:,3))/2;
% ave2 = (crs_subfram_res(:,2) + crs_subfram_res(:,4))/2;
% ave_ = [transpose(ave1);transpose(ave2)];
% ave_inter = reshape(ave_,400,1);
% % Frequency averaging
% window_size = 0;
% fre_aver = zeros(400,1);
% for n = 1: 400
%     if n < 10
%         fre_aver(n) = mean(ave_inter(n-window_size:n+window_size));
%         window_size = window_size +1;
%     elseif 10 <= n && n <= 391
%         fre_aver(n) = mean(ave_inter(n-window_size:n+window_size));
%     else
%         window_size = window_size - 1;
%         fre_aver(n) = mean(ave_inter(n-window_size:n+window_size));
%     end
% end

% 1200-> 600 199->99 1194->594 200->100 796->396 1203->603 1197->597
% frequency Interpolation
grid_estimation = zeros(600,14);
quotient_c = zeros(99,4);
h_est_coding_reshape = reshape(h_est_coding,100,4);
h_est_coding_reshape(:,1) = h_est_coding_reshape(:,1);

h_est_coding_reshape_22 = reshape(h_est_coding_reshape(1:99,:),1,396);

for n =1:99
    quotient_c(n,:) = (h_est_coding_reshape(n+1,:) - h_est_coding_reshape(n,:))/6;
end
quotient_c_reshape = reshape(quotient_c,1,396);

h_est_coding_reshape_add_one = zeros(6,396);
for n = 1:6
    h_est_coding_reshape_add_one(n,:) = h_est_coding_reshape_22 + (n-1) * quotient_c_reshape;
end
h_est_coding_reshape_add_one_1 = reshape(h_est_coding_reshape_add_one,594,4);


h_est_insert_re_Re_ = zeros(603,4);
h_est_coding_reshape_insert_re_Re_ =  [h_est_coding_reshape_add_one_1;repmat(h_est_coding_reshape_add_one_1(end,:),6,1)];
h_est_insert_re_Re_(4:end,:) = h_est_coding_reshape_insert_re_Re_;

for n = 1:6
    h_est_insert_re_Re_(597+n,:) = h_est_coding_reshape_insert_re_Re_(594+n,:) + n * quotient_c(99,:);
end

h_est_insert_re_Re_(1:3,:) = repmat(h_est_insert_re_Re_(4,:),3,1);

h_est_insert_re_Re_(3,:) = h_est_insert_re_Re_(3,:) - 1 * quotient_c(1,:);
h_est_insert_re_Re_(2,:) = h_est_insert_re_Re_(2,:) - 2 * quotient_c(1,:);
h_est_insert_re_Re_(1,:) = h_est_insert_re_Re_(1,:) - 3 * quotient_c(1,:);


%% time interpolation
grid_estimation(:,1:4) = repmat(h_est_insert_re_Re_(4:end,1),1,4);
grid_estimation(:,5:7) = repmat(h_est_insert_re_Re_(1:600,2),1,3);
grid_estimation(:,8:11) = repmat(h_est_insert_re_Re_(4:end,3),1,4);
grid_estimation(:,12:14) = repmat(h_est_insert_re_Re_(1:600,4),1,3);

quotient_2 = zeros(600,3);
quotient_2(:,1) = (grid_estimation(:,5) - grid_estimation(:,1))/4;
quotient_2(:,2) = (grid_estimation(:,8) - grid_estimation(:,5))/3;
quotient_2(:,3) = (grid_estimation(:,12) - grid_estimation(:,8))/4;

grid_estimation_cc = zeros(600,14);
for n = 1:3
    grid_estimation(:,n+1) = grid_estimation(:,n+1) + n * quotient_2(:,1);
end
for n = 6:7
    grid_estimation(:,n) = grid_estimation(:,n) + (n-5) * quotient_2(:,2);
end
for n = 9:11
    grid_estimation(:,n) = grid_estimation(:,n) + (n-8) * quotient_2(:,3);
end
for n = 13:14
    grid_estimation(:,n) = grid_estimation(:,n) + (n-12) * quotient_2(:,3);
end

figure(31)
surf(abs(grid_estimation(:,:,1,1)/32))  %%channel_status_grid  grid_estimation
% view(0,0)
% figure(32)
% surf(abs(channel_status_grid(:,:,1,1)))  %%channel_status_grid  grid_estimation
%
% difference = grid_estimation/32 - channel_status_grid;
% figure(38)
% plot(abs(difference(:,7)));

%% Channel Compensation
% eqGrid = lteEqualizeMMSE(rxGrid, estChannel, noiseEst);
% [eqGrid,csi] = lteEqualizeZF(rxGrid,estChannel);

eqGrid_z = rxGrid./grid_estimation;

eqGrid_z_re = reshape(eqGrid_z,numel(eqGrid_z),1);
real_grid =real(eqGrid_z_re);
imag_grid =imag(eqGrid_z_re);
figure(22)
% hold on
scatter(real_grid,imag_grid)

%% Channel Estimator Configuration
cec.PilotAverage = 'UserDefined'; % Pilot averaging method
cec.FreqWindow = 9;               % Frequency averaging window in REs
cec.TimeWindow = 9;               % Time averaging window in REs
cec.InterpType = 'Cubic';         % Cubic interpolation
cec.InterpWinSize = 3;            % Interpolate up to 3 subframes
% simultaneously
cec.InterpWindow = 'Centred';     % Interpolation windowing method

%% Channel Estimation
% rxGrid(cellRsInd) = test;
enb.NSubframe = 0;
[estChannel, noiseEst] = lteDLChannelEstimate(enb,cec,rxGrid); %%fftc_output_grid  rxGrid
% rxGrid()
figure()
scatter(real(rxGrid(1:6:end,1)), imag(rxGrid(1:6:end,1)));
hold on;
scatter(real(rxGrid(4:6:end,1)), imag(rxGrid(4:6:end,1)));
figure(32)
surf(abs(estChannel(:,:,1,1)))  %%channel_status_grid estChannel

%% MMSE Equalization
eqGrid = lteEqualizeZF(rxGrid, estChannel);
scatterplot(eqGrid(1:6:end,1));
% scatterplot(eqGrid(1:12:end,1));
combbb = reshape(eqGrid,numel(eqGrid),1);

real_ =real(combbb);
imag_ =imag(combbb);
figure()
% hold on
scatter(real_,imag_)

% figure()
% scatter(real_(601:1200),imag_(601:1200))
% figure()
% scatter(real_(1201:1800),imag_(1201:1800))
% figure()
% scatter(real_(1801:2400),imag_(1801:2400))