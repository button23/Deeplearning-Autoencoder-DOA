close all;
clear;
clc;

mcss = [0, 1, 3, 5, 8];    % BPSK:0  QPSK:1  16QAM:3  64QAM:5  256QAM: 8
dataLength = [1, 2, 4, 8, 12]; % BPSK:1  QPSK:2  16QAM:4  64QAM:8  256QAM:12
numTrain = 1000;
numTest = 100;
validLength = 80;


snrr_record = zeros(numTest,length(mcss));
train_data_before = zeros(numTrain * length(mcss),80);
train_data = zeros(numTrain * length(mcss),52);
test_data_before = zeros(numTest * length(mcss),80);
test_data = zeros(numTest * length(mcss),52);


for i=1:length(mcss)
    %% Create a phase and frequency offset object and introduce a 2 Hz frequency offset.
    %% wlan waveform
    cfgVHT = wlanVHTConfig('ChannelBandwidth','CBW20','MCS',mcss(i));

    cfgVHT.APEPLength = 50 * dataLength(i);
    
    % Random snr
    % snrr = 5;
    snrr = randi([-5,20], (numTrain+numTest),1);
    % Record the SNR of the test data. 
    snrr_record(:,i) = snrr((numTrain+1):end);
    
    tgacChan = wlanTGacChannel('SampleRate',20e6,'ChannelBandwidth','CBW20','DelayProfile','Model-B');

    % pfOffset = comm.PhaseFrequencyOffset('SampleRate',20e6,'FrequencyOffsetSource','Input port');
    
    %%
    for mjk =1:(numTrain + numTest)
        %% Data Convertion
        % Create random bit stream
        %       rng(25)
        PSDU = randi([0 1],cfgVHT.PSDULength*8,1);
        wholeWaveform = wlanWaveformGenerator(PSDU,cfgVHT);
        
        % Channel and AWGN
        chNoise = comm.AWGNChannel('NoiseMethod','Signal to noise ratio (SNR)',...
        'SNR', snrr(i));
        rxSig = chNoise(tgacChan(wholeWaveform));  %% wholeWaveform_scaling  wholeWaveform
%         rxSig = wholeWaveform;  %% wholeWaveform_scaling  wholeWaveform
        
        %% Synchronization
        [startOffset,meric] = wlanSymbolTimingEstimate(rxSig,'CBW20'); %% rxSig  RX_Data WIFI_Waveform
        RX_Data_ = rxSig(startOffset+1:end);
        if mjk <= numTrain
          train_data_before(((i-1)*numTrain+mjk),:) = RX_Data_(801:(800+validLength));
        else
          test_data_before(((i-1)*numTest+mjk-numTrain),:) = RX_Data_(801:(800+validLength));
        end
    end
end

%% Training set processing
for ii = 1:numTrain*length(mcss)
    symbb=train_data_before(ii,:);
    data_d_s = symbb(17:end);
    dsp_test_data_before = fft(data_d_s);
    dsp_test_fftshift = fftshift(dsp_test_data_before,1) * sqrt(56) * 1/64;
    zeropadding_remove_origin = [dsp_test_fftshift(5:32),dsp_test_fftshift(34:61)];
    pilot_remove =[zeropadding_remove_origin(1:7),zeropadding_remove_origin(9:21),zeropadding_remove_origin(23:34),zeropadding_remove_origin(36:48),zeropadding_remove_origin(50:56)];
    train_data(ii,: ) =  pilot_remove.';
end

%% Test set processing
for pp = 1:numTest*length(mcss)
    symbb=test_data_before(pp,:);
    data_d_s = symbb(17:end);
    dsp_test_data_before = fft(data_d_s);
    dsp_test_fftshift = fftshift(dsp_test_data_before,1) * sqrt(56) * 1/64;
    zeropadding_remove_origin = [dsp_test_fftshift(5:32),dsp_test_fftshift(34:61)];
    pilot_remove =[zeropadding_remove_origin(1:7),zeropadding_remove_origin(9:21),zeropadding_remove_origin(23:34),zeropadding_remove_origin(36:48),zeropadding_remove_origin(50:56)];
    test_data(pp,: ) =  pilot_remove.';
end

snrr_record_reshape = reshape(snrr_record,numel(snrr_record),1);
test_snr = snrr_record_reshape;

% Training set label
[hight_d, length_d] = size(train_data_before);
train_label = zeros(hight_d,1);
for ii = 1 : length(mcss)
    train_label((hight_d/length(mcss) * (ii-1)+1):hight_d/length(mcss) * ii) = ii-1;
end

% Test set label
[hight_d, length_d] = size(test_data_before);
test_label = zeros(hight_d,1);
for ii = 1 : length(mcss)
    test_label((hight_d/length(mcss) * (ii-1)+1):hight_d/length(mcss) * ii) = ii-1;
end

% %
% train_data_before___ = train_data_before__(:,1:1000);
% %%
% anss = train_data_before(:,10);
% data_ = anss(801:end);
% data_d =reshape(data_,80,4);
% data_d_s = data_d(17:end,:);
%
% dsp_test_data_before = fft(data_d_s);
% dsp_test_fftshift = fftshift(dsp_test_data_before,1) * sqrt(56) * 1/64;
% dsp_test_fftshift_ = reshape(dsp_test_fftshift,numel(dsp_test_fftshift),1);
%
% figure(6)
% xlim([-1.5 1.5])
% ylim([-1.5 1.5])
% pilot_remove_ =(reshape(dsp_test_fftshift_,[64*4,1]));
% sz = 25;
% scatter(real(pilot_remove_),imag(pilot_remove_),sz,'*')
% xlabel('In-phase')
% ylabel('Quadrature')
% title('Received 256 QAM Constellation')
%
% %% CONCATE
% train_data_before = zeros(10000*5,1000);
% train_data_before(1:10000,:) = BPSK(1:1000,:).';
% train_data_before(10001:20000,:) = QPSK(1:1000,:).';
% train_data_before(20001:30000,:) = QAM16(1:1000,:).';
% train_data_before(30001:40000,:) = QAM64(1:1000,:).';
% train_data_before(40001:50000,:) = QAM256(1:1000,:).';
%

% label_data(20001:300) = 2;
% label_data(30001:400) = 3;
% label_data(40001:500) = 4;





