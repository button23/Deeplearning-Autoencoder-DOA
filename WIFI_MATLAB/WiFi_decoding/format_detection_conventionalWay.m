clc
% clear

%% Original way
sr = 20e6;
Tlstf = 8e-6;
Tlltf = 8e-6;
cfgVHT = wlanVHTConfig('ChannelBandwidth','CBW20','MCS',0); % bpsk 1/2 data rate 3

idxlltf = Tlstf*sr+(1:Tlltf*sr);
TdetectionSymbols = 12e-6;
idxDetectionSymbols = (Tlstf+Tlltf)*sr+(1:TdetectionSymbols*sr);

index_discard = zeros(length(test_label),1);
kk1 = 1;
kk2 = 1;
% format_record = cell(length(test_data(:,1)),1);
format_record = zeros(length(test_data(:,1)),1);


for i = 1:length(test_data(:,1))
    dataTT = test_data(i,:);
    %Syncronization
    startOffset = wlanSymbolTimingEstimate(dataTT.', ...
        cfgVHT.ChannelBandwidth);
    if (isempty(startOffset))
        index_discard(kk1) = i;
        kk1 = kk1+1;
        continue
    end
    if startOffset< 0 || startOffset> 80
        index_discard(kk1) = i;
        kk1 = kk1+1;
        continue
    end
    
    rxSig_sync = dataTT(startOffset+1:startOffset+560).';
    
    oneFrame_ltf = rxSig_sync(idxlltf,:);
    
    lltfDemod = wlanLLTFDemodulate(oneFrame_ltf,'CBW20');
    
    chEst = wlanLLTFChannelEstimate(lltfDemod,'CBW20');
    
    noiseVarEst = 10.^(-test_snr(i)/20);
    
    
    in = rxSig_sync(idxDetectionSymbols,:);
    
    cfgRec = wlanRecoveryConfig('OFDMSymbolOffset',0.5,...
        'PilotPhaseTracking','None');
    
    format_detect = wlanFormatDetect(in,chEst,noiseVarEst,'CBW20',cfgRec);
    
    %     format_record{i} = format_detect;
    if strcmp(format_detect, 'VHT')
        format_record(kk2) = 0;
    elseif strcmp(format_detect,  'Non-HT')
        format_record(kk2) = 2;
    else
        format_record(kk2) = 1;
    end
    kk2 = kk2+1;
end
format_record(kk2:end) = [];
test_label(index_discard(1:kk1-1)) = [];
test_snr(index_discard(1:kk1-1)) = [];

%%
SNRR = -10:20;
snr_cc = zeros(31,1);
snr_dis_ = zeros(31,1);

erra = test_label - format_record;

errindd = find(erra~=0);
error_snr = test_snr(errindd);

for pp = 1:length(SNRR)
    snr_cc(pp) = length(find(error_snr == SNRR(pp)));
end
for ppp = 1:length(SNRR)
    snr_dis_(ppp) = length(find(test_snr == SNRR(ppp)));
end
snr_dis = 3000 - snr_dis_;
snrfinal = snr_dis+snr_cc;

figure(1)
plot(SNRR,(1-snrfinal/3000));

