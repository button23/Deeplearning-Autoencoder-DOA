%% Simulate Initialization
clc
clear; close all;
% the reason why control channels including PBCH, PCFICH, PDCCH, etc are
% divided with sqrt(2) is to normalize Tx antenna powers
% based on RMC waveform R.10
% decoding procedures for control channels should be implemented for future

%% LTE Transmition Mode 3(MIMO, 2 by 2, 10MHz, 64QAM) eNB Setting
txMode = 'TM3';
MCS = '64QAM';

TM3_Setting = []; % clear simulationParameters
TM3_Setting.NDLRB = 50;
TM3_Setting.PDSCH.TargetCodeRate = 0.75;
TM3_Setting.PDSCH.PRBSet = (0:TM3_Setting.NDLRB-1)';
TM3_Setting.CFI = 1;
TM3_Setting.TotSubframes = 10; % Generate one subframe at a time
TM3_Setting.NSubframe = 0; % 0 ~9

% TM3_Setting.PDSCH.TxScheme = 'CDD';
TM3_Setting.PDSCH.TxScheme = 'SpatialMux';
TM3_Setting.PDSCH.DCIFormat = 'Format2A';
TM3_Setting.CellRefP = 2;
TM3_Setting.PDSCH.Modulation = {MCS, MCS};
TM3_Setting.PDSCH.NLayers = 2;

nCodewords = 2;
eNB = lteRMCDL(TM3_Setting);
eNB.PDSCH.RVSeq = [0,1,2,3];

%% LTE Transmition Mode 3(MIMO, 2 by 2, 10MHz, 64QAM) Parameter Setting
Num_Subframe    = eNB.NSubframe;
numREperOFDMSym   = (eNB.NDLRB) * 12;
numREperSubframe = (numREperOFDMSym * 14) * (eNB.CellRefP);

Tx_Antenna           = eNB.CellRefP;
TransBlock_Size      = eNB.PDSCH.TrBlkSizes(1,Num_Subframe+1);
CodedTransBlock_Size = eNB.PDSCH.CodedTrBlkSizes(1,Num_Subframe+1);
Scramble_Sequence    = ltePDSCHPRBS(eNB,eNB.PDSCH.RNTI,0,CodedTransBlock_Size);
Descramble_Sequemce  = 1-2*Scramble_Sequence; % 0 ¡ú 1, 1 ¡ú -1

%% LTE Transmition Mode 3(MIMO, 2 by 2, 10MHz, 64QAM) PDSCH (Physical Downlink Shared Channel) Encoding
Test_Num = 1;

Information_bit = zeros(TransBlock_Size, Tx_Antenna);
switch (Test_Num)
    case 1
        % DSP, RMC Waveform Pattern Bits        
        Test_bit = [1 1 0 0]; Test_bit2 = [1 1 0 0];
        for idx = 1:1:TransBlock_Size/4 
            for i = 1:1:4
                Information_bit(i+4*(idx-1),1) = Test_bit(i);
                Information_bit(i+4*(idx-1),2) = Test_bit2(i);
            end
        end
        
    case 2
        % Random Bits
        Information_bit = randi([0 1], TransBlock_Size, 2);        
end

RateMatch_Data = zeros(CodedTransBlock_Size,   Tx_Antenna);
Scramble_Data  = zeros(CodedTransBlock_Size,   Tx_Antenna);
PDSCH_Data     = zeros(CodedTransBlock_Size/6, Tx_Antenna);

for Ant = 1:1:Tx_Antenna
    
    % CRC insertion(Matlab Function Call)
    CRC_insertion_Data = lteCRCEncode(Information_bit(:,Ant),'24A');
    
    % Segmentation(Matlab Function Call)
    Segment_Data = lteCodeBlockSegment(CRC_insertion_Data);
    
    % Turbo Encode(Matlab Function Call)
    Turbo_Encode_Data = lteTurboEncode(Segment_Data);
    
    % RateMatch(Matlab Function Call)
    RateMatch_Data(:,Ant) = lteRateMatchTurbo(Turbo_Encode_Data, CodedTransBlock_Size, 0, eNB.PDSCH);
    
    % Scramble
    Scramble_Data(:,Ant) = xor(Scramble_Sequence,RateMatch_Data(:,Ant));
    
    % 64QAM Modulation
    PDSCH_Data(:,Ant) = lteSymbolModulate(Scramble_Data(:,Ant),MCS);
end


%% LTE Transmition Mode 3(MIMO, 2 by 2, 10MHz, 64QAM) - Precoding
Matalb_Precoded = lteDLPrecode(eNB, eNB.PDSCH, PDSCH_Data);
Matalb_Precoded = round(Matalb_Precoded,4);

Precod_PDSCH = zeros(length(PDSCH_Data),Tx_Antenna);
if (eNB.PDSCH.TxScheme == "CDD")
    for i= 1:1:length(PDSCH_Data)
        Precod_PDSCH(i,1) = (PDSCH_Data(i,1)+PDSCH_Data(i,2))/2;
    end
    
    for i = 1:1:length(PDSCH_Data)/2
        Precod_PDSCH(2*i-1,2) = (PDSCH_Data(2*i-1,1)-PDSCH_Data(2*i-1,2))/2;
        Precod_PDSCH(2*i,2)   = (-PDSCH_Data(2*i,1)+PDSCH_Data(2*i,2))/2;
    end
    
elseif(eNB.PDSCH.TxScheme == "SpatialMux") % use codebook index 1 => 0.5*[1 1;1 -1] default codebook index is 1 for your information
    for i= 1:1:length(PDSCH_Data)
        Precod_PDSCH(i,1) = (PDSCH_Data(i,1) + PDSCH_Data(i,2))/2;
        Precod_PDSCH(i,2) = (PDSCH_Data(i,1) - PDSCH_Data(i,2))/2;
    end    
end

Precod_PDSCH   = round(Precod_PDSCH,4);
Error_Precode_PDSCH = sum(Matalb_Precoded ~= Precod_PDSCH);

%% Check Matalb Function DLSCH, PDSCH
[Matlab_DLSCH1, Matalb_DLSCH1_Info] = lteDLSCH(eNB, eNB.PDSCH, eNB.PDSCH.CodedTrBlkSizes(1), Information_bit(:,1));
[Matlab_DLSCH2, Matalb_DLSCH2_Info] = lteDLSCH(eNB, eNB.PDSCH, eNB.PDSCH.CodedTrBlkSizes(1), Information_bit(:,2));

Error_DLSCH(1,1) = sum(RateMatch_Data(:,1) ~= Matlab_DLSCH1);
Error_DLSCH(1,2) = sum(RateMatch_Data(:,2) ~= Matlab_DLSCH2);

% %% DSP Encoder Output Check
% DSP_EncoderScaling = 2^10; % for BCP input/output scaling
% 
% DSP_Precode_PDSCH = DSP_Modulate(Precod_PDSCH, DSP_EncoderScaling);

%% LTE Transmition Mode 2(Tx Diversity, 2 by 1, 10MHz, 64QAM) Reference Signal
CellRS     = lteCellRS(eNB);
CellRS_idx = lteCellRSIndices(eNB);

LTE_Symbol = zeros(8400*Tx_Antenna, 1);
for idx = 1:1:length(CellRS_idx)
    LTE_Symbol(CellRS_idx(idx)) = CellRS(idx);
end

CellRS_idx2 = zeros(length(CellRS_idx)/2,1);
for i = 1:1:length(CellRS_idx)/8 % please refer to LTE grid when there are two antenna ports. +-3 means CRS offset of different antenna ports
    CellRS_idx2(i) = CellRS_idx(i)+3;
    CellRS_idx2(i+100) = CellRS_idx(i+100)-3;
    CellRS_idx2(i+200) = CellRS_idx(i+200)+3;
    CellRS_idx2(i+300) = CellRS_idx(i+300)-3;
end

%% LTE Transmition Mode 2(Tx Diversity, 2 by 1, 10MHz, 64QAM) Sync PSS, SSS

if (Num_Subframe == 0 || Num_Subframe == 5)
    PSS = ltePSS(eNB); SSS = lteSSS(eNB);
    PSS_idx = zeros(length(PSS), 1); SSS_idx = zeros(length(SSS), 1);
    
    for idx = 1:1:1
        PSS_idx(:, idx) = ltePSSIndices(eNB,idx-1);
        SSS_idx(:, idx) = lteSSSIndices(eNB,idx-1);
    end
    
    for Sub = 1:1:1
        for idx = 1:1:length(PSS_idx)
            LTE_Symbol(PSS_idx(idx, Sub)) = PSS(idx);
            LTE_Symbol(SSS_idx(idx, Sub)) = SSS(idx);            
        end
    end
end

%% LTE Transmition Mode 3(MIMO, 2 by 2, 10MHz, 64QAM) Control Channl - PBCH
if (Num_Subframe == 0)
    % PBCH Information
%     MIB         = flipud(lteMIB(eNB));
    MIB      = lteMIB(eNB);
    MIB_CRC  = lteCRCEncode(MIB,'16');
    MIB_Conv = lteConvolutionalEncode(MIB_CRC);
    MIB_RM   = lteRateMatchConvolutional(MIB_Conv,1920);
    PBCH_Scramble_Sequence = ltePBCHPRBS(eNB,1920);
    MIB_Scramble = xor(MIB_RM,PBCH_Scramble_Sequence);
    PBCH_Symbol = lteSymbolModulate(MIB_Scramble, 'QPSK');
    PBCH_idx    = ltePBCHIndices(eNB, {'1based'});
    
    Precod_PBCH = zeros(length(PBCH_Symbol),Tx_Antenna);
    if (eNB.PDSCH.TxScheme == "CDD")
        for idx = 1:1:length(PBCH_Symbol)
            Precod_PBCH(idx,1) = PBCH_Symbol(idx)/sqrt(2); 
        end
        
        for idx = 1:1:length(PBCH_Symbol)/2
            Precod_PBCH(2*idx-1,2) = -1*conj(Precod_PBCH(2*idx,1));
            Precod_PBCH(2*idx,2) = conj(Precod_PBCH(2*idx-1,1));
        end
    elseif(eNB.PDSCH.TxScheme == "SpatialMux")
        for i= 1:1:length(PBCH_Symbol)
            Precod_PBCH(i,1) = PBCH_Symbol(i,1)/sqrt(2);
            Precod_PBCH(i,2) = PBCH_Symbol(i,1)/sqrt(2);
        end
    end    
    
    for Sub = 1:1:Tx_Antenna
        for idx = 1:1:length(PBCH_idx)
            LTE_Symbol(PBCH_idx(idx, Sub), 1) = Precod_PBCH(idx, Sub);
        end
    end
end

Matalb_Precoded_PBCH = lteDLPrecode(eNB, eNB.PDSCH, PBCH_Symbol);
Error_Precode_PBCH = sum(round(Precod_PBCH,4) ~= round(Matalb_Precoded_PBCH,4));

%% LTE Transmition Mode 3(MIMO, 2 by 2, 10MHz, 64QAM) Control Channl - PCFICH
if (Num_Subframe == 0)
    PCFICH_CW     = lteCFI(eNB);
    PCFICH_Symbol = ltePCFICH(eNB, PCFICH_CW);
    PCFICH_idx    = ltePCFICHIndices(eNB, {'1based'});
    
    for Sub = 1:1:Tx_Antenna
        for idx = 1:1:length(PCFICH_idx)
            LTE_Symbol(PCFICH_idx(idx, Sub), 1) = PCFICH_Symbol(idx, Sub);
        end
    end
end

%% LTE Transmition Mode 3(MIMO, 2 by 2, 10MHz, 64QAM) Control Channl - PCFICH Precode
if (Num_Subframe == 0)
    PCFICH_Scramble_Sequence = ltePCFICHPRBS(eNB, length(PCFICH_CW));
    PCFICH_Scramble = xor(PCFICH_CW,PCFICH_Scramble_Sequence);
    PCFICH_Data = lteSymbolModulate(PCFICH_Scramble, 'QPSK');
    
    Precod_PCFICH = zeros(length(PCFICH_Data),Tx_Antenna);
    if (eNB.PDSCH.TxScheme == "CDD")
        for idx = 1:1:length(PCFICH_Data)
            Precod_PCFICH(idx,1) = PCFICH_Data(idx)/sqrt(2);
        end
        
        for idx = 1:1:length(PCFICH_Data)/2
            Precod_PCFICH(2*idx-1,2) = -1*conj(Precod_PCFICH(2*idx,1));
            Precod_PCFICH(2*idx,2) = conj(Precod_PCFICH(2*idx-1,1));
        end
    elseif(eNB.PDSCH.TxScheme == "SpatialMux")
        for i= 1:1:length(PCFICH_Data)
            Precod_PCFICH(i,1) = (PCFICH_Data(i,1))/sqrt(2);
            Precod_PCFICH(i,2) = (PCFICH_Data(i,1))/sqrt(2);
        end
    end
end

Matalb_Precoded_PCFICH = lteDLPrecode(eNB, eNB.PDSCH, PCFICH_Data);
Error_Precode_PCFICH = sum(round(Precod_PCFICH,4) ~= round(Matalb_Precoded_PCFICH,4));
 
%% LTE Transmition Mode 2(Tx Diversity, 2 by 2, 10MHz, 64QAM) Control Channl - PHICH
if (Num_Subframe == 0)    
    [PHICH_Symbol,PHICH_Info] = ltePHICH(eNB,eNB.HISet);
    PHICH_idx = ltePHICHIndices(eNB, {'1based'});
    
    for Sub = 1:1:Tx_Antenna
        for idx = 1:1:length(PHICH_idx)
            LTE_Symbol(PHICH_idx(idx,Sub), 1) = PHICH_Symbol(idx, Sub);
        end
    end    
end

%% LTE Transmition Mode 2(Tx Diversity, 2 by 2, 10MHz, 64QAM) Control Channl - PDCCH
DCI = lteDCI(eNB, eNB.PDSCH);
DCI_Info = lteDCIInfo(eNB, eNB.PDSCH);
DCI_MessageBits = [0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 1 0 0 0 1 0 0 0 0];
DCI_MessageBits_length = DCI_Info.Format1;

DCI_CRC  = lteCRCEncode(DCI_MessageBits,'16', eNB.PDSCH.RNTI);
DCI_Conv = lteConvolutionalEncode(DCI_CRC);
DCI_Rm   = lteRateMatchConvolutional(DCI_Conv, 288);

% lteDCIEncode
DCI_Codeword = lteDCIEncode(eNB.PDSCH, DCI_MessageBits);
Error_DCI = sum(DCI_Codeword ~= DCI_Rm);

Candidate    = ltePDCCHSpace(eNB, eNB.PDSCH);
PDCCH_Info   = ltePDCCHInfo(eNB);
PDCCH_Symbol = zeros(PDCCH_Info.MTot/2,1);

PDCCH_Scramble_Sequence    = ltePDCCHPRBS(eNB,PDCCH_Info.MTot);
PDCCH_Descramble_Sequemce  = 1 - 2*PDCCH_Scramble_Sequence(Candidate(2,1):Candidate(2,2)); % 0 ¡ú 1, 1 ¡ú -1

DCI_Scramble = xor(PDCCH_Scramble_Sequence(Candidate(2,1):Candidate(2,2)), DCI_Rm);
DCI_Modulate = lteSymbolModulate(DCI_Scramble,'QPSK');

PDCCH_Symbol(1:144) = DCI_Modulate;
PDCCH_Interleave = ltePDCCHInterleave(eNB,PDCCH_Symbol);
    
Precod_PDCCH = zeros(length(PDCCH_Interleave),Tx_Antenna);
if (eNB.PDSCH.TxScheme == "CDD")
    for idx = 1:1:length(PDCCH_Interleave)
        Precod_PDCCH(idx,1) = PDCCH_Interleave(idx)/sqrt(2);
    end
    
    for idx = 1:1:length(PDCCH_Interleave)/2
        Precod_PDCCH(2*idx-1,2) = -1*conj(Precod_PDCCH(2*idx,1));
        Precod_PDCCH(2*idx,2) = conj(Precod_PDCCH(2*idx-1,1));
    end
elseif(eNB.PDSCH.TxScheme == "SpatialMux")
    for i= 1:1:length(PDCCH_Interleave)
        Precod_PDCCH(i,1) = (PDCCH_Interleave(i,1))/sqrt(2);
        Precod_PDCCH(i,2) = (PDCCH_Interleave(i,1))/sqrt(2);
    end
end

Matalb_Precoded_PDCCH = lteDLPrecode(eNB, eNB.PDSCH, PDCCH_Interleave);
Error_Precode_PDCCH = sum(round(Precod_PDCCH,4) ~= round(Matalb_Precoded_PDCCH,4));

%% LTE Transmition Mode 3(MIMO, 2 by 2, 10MHz, 64QAM) Control Channl - Layer Mapping
PDSCH_idx = ltePDSCHIndices(eNB,eNB.PDSCH, ...
    eNB.PDSCH.PRBSet,{'1based','ind'});

for Sub = 1:1:Tx_Antenna
    for idx = 1:1:length(Precod_PDSCH)
        LTE_Symbol(PDSCH_idx(idx,Sub), 1) = Precod_PDSCH(idx, Sub);
    end
end

% for i = 1 : 10
LTE_Mapping = zeros(600,14,Tx_Antenna);
for Ant = 1:1:Tx_Antenna
    for Sub = 1:1:14
        for idx = 1:1:600
            LTE_Mapping(idx, Sub, Ant) = LTE_Symbol(idx+600*(Sub-1) + 8400*(Ant-1));
        end
    end
end

Antenna_Mapping1 = LTE_Mapping(:,:,1);
Antenna_Mapping2 = LTE_Mapping(:,:,2);

% TDD (#1 #2 subframe ->0)

%% LTE Transmition Mode 3(MIMO, 2 by 2, 10MHz, 64QAM) - IFFT, CP
    [LTE_Waveform, LTE_Waveform_Info] = lteOFDMModulate(eNB, LTE_Mapping);
%   if i == 2 || i == 3   
%     % fdd
%     LTE_Waveform_tdd((i-1)*15360+1:15360*i,:) = 0;
%   else
%           LTE_Waveform_tdd((i-1)*15360+1:15360*i,:) = LTE_Waveform;
%   end
  
%  end
%% LTE Transmition Mode 3(MIMO, 2 by 2, 10MHz, 64QAM) - Fading Channel

%% [DSP] LTE Transmition Mode 3(MIMO, 2 by 2, 10MHz, 64QAM) - Control Channel
DSP_EncoderScaling = 2^10;

if (Num_Subframe == 0)
    DSP_PBCH = DSP_Modulate(Precod_PBCH, DSP_EncoderScaling);
    DSP_PHICH = DSP_Modulate(PHICH_Symbol, DSP_EncoderScaling);
    DSP_PDCCH = DSP_Modulate(Precod_PDCCH, DSP_EncoderScaling);
    DSP_PCFICH = DSP_Modulate(PCFICH_Symbol, DSP_EncoderScaling);
end

%% [DSP] LTE Transmition Mode 3(MIMO, 2 by 2, 10MHz, 64QAM) - GRC Sync, DSP Input Setting
DSP_OutputScaling = 2^18;

LTE_Waveform_Check(1:15360,1) = LTE_Waveform(:,1);
LTE_Waveform_Check(15361:30720,1) = LTE_Waveform(:,2);

DSP_Antenna1 = DSP_Modulate(LTE_Waveform(:,1), DSP_OutputScaling);
DSP_Antenna2 = DSP_Modulate(LTE_Waveform(:,2), DSP_OutputScaling);

DSP_LTE_Waveform = DSP_Modulate(LTE_Waveform, DSP_OutputScaling);

for sub = 1:1:30
    for idx = 1:1:512
        DSP_MatlabInput(sub,idx)  = DSP_LTE_Waveform(idx+512*(sub-1),1);
        DSP_MatlabInput2(sub,idx) = DSP_LTE_Waveform(idx+512*(sub-1),2);
    end
end

%% [DSP] LTE Transmition Mode 3(MIMO, 2 by 2, 10MHz, 64QAM) - DSP Input Setting
DSP_TestOn = 0; % 0 : No Fading Channel // 1 : RF Test // 2 : Fading Channel

if (DSP_TestOn == 0)
    
    corrcfg_PSS.PSS = 'On';    
    [offset_PSS,corr_PSS] = lteDLFrameOffset(eNB,LTE_Waveform,corrcfg_PSS);
    
%     figure(1)
%     plot(abs(corr_PSS))
end

if(DSP_TestOn == 1)
    for idx = 1:1:length(VarName1)/2
        DSP_signal1(idx,1) = VarName1(2*idx) + 1j*(VarName1(2*idx-1));
        DSP_signal2(idx,1) = VarName2(2*idx) + 1j*(VarName2(2*idx-1));
    end    
    
    corrcfg.PSS = 'On';
    [offset1,corr1] = lteDLFrameOffset(eNB,DSP_signal1,corrcfg);
    [offset2,corr2] = lteDLFrameOffset(eNB,DSP_signal2,corrcfg);
    
%     figure(1)
%     plot(abs(corr1))
%     
%     figure(2)
%     plot(abs(corr2))
    
    RF_Subframe1_Data1 = DSP_signal1((offset1+1) : (offset1+15360));
    RF_Subframe1_Data2 = DSP_signal2((offset1+1) : (offset1+15360));
     
    for sub = 1:1:30
        for idx = 1:1:512
            DSP_RFInput(sub,idx)  = RF_Subframe1_Data1(idx+512*(sub-1));
            DSP_RFInput2(sub,idx) = RF_Subframe1_Data2(idx+512*(sub-1));
        end
    end
end

%% LTE Transmition Mode 3(MIMO, 2 by 2, 10MHz, 64QAM) - Recive Antenna 1
switch (DSP_TestOn)
    case 0
        Recive_Mappings_Ant1 = lteOFDMDemodulate(eNB, LTE_Waveform);
        Recive_Mapping_Ant1  = Recive_Mappings_Ant1(:, :, 1);
    case 1
        Recive_Mapping_Ant1 = lteOFDMDemodulate(eNB, RF_Subframe1_Data1);
end

Recive_LTE_Ant1 = zeros(numREperSubframe/2,1);
for Sub = 1:1:14
    for idx = 1:1:numREperOFDMSym
        Recive_LTE_Ant1(idx + 600*(Sub-1)) = Recive_Mapping_Ant1(idx, Sub);
    end
end

%% LTE Transmition Mode 3(MIMO, 2 by 2, 10MHz, 64QAM) - Recive Antenna 2
switch (DSP_TestOn)
    case 0
        Recive_Mappings_Ant2 = lteOFDMDemodulate(eNB, LTE_Waveform);
        Recive_Mapping_Ant2  = Recive_Mappings_Ant2(:, :, 2);
    case 1
        Recive_Mapping_Ant2 = lteOFDMDemodulate(eNB, RF_Subframe1_Data2);
end

Recive_LTE_Ant2 = zeros(numREperSubframe/2,1);
for Sub = 1:1:14
    for idx = 1:1:numREperOFDMSym
        Recive_LTE_Ant2(idx + 600*(Sub-1)) = Recive_Mapping_Ant2(idx, Sub);
    end
end

%% Recive Signal Channel Reference Estimation 
 
Reference_Estimate1 = zeros(length(CellRS_idx)/2,1);
Reference_Estimate2 = zeros(length(CellRS_idx)/2,1);
for idx = 1:1:length(CellRS_idx)/2
    Reference_Estimate1(idx) = Recive_LTE_Ant1(CellRS_idx(idx))/CellRS(idx);
    Reference_Estimate2(idx) = Recive_LTE_Ant1(CellRS_idx2(idx))/CellRS(idx);
end

Reference_Estimate3 = zeros(length(CellRS_idx)/2,1);
Reference_Estimate4 = zeros(length(CellRS_idx)/2,1);
for idx = 1:1:length(CellRS_idx)/2
    Reference_Estimate3(idx) = Recive_LTE_Ant2(CellRS_idx(idx))/CellRS(idx);
    Reference_Estimate4(idx) = Recive_LTE_Ant2(CellRS_idx2(idx))/CellRS(idx);
end

%% Recive Signal Channel Estimation 

[Estimate_Channel1,Estimates1_grid] = Channel_Estimate_odd(Reference_Estimate1,numREperSubframe/2, CellRS_idx);
[Estimate_Channel2,Estimates2_grid] = Channel_Estimate_even(Reference_Estimate2,numREperSubframe/2, CellRS_idx2);
[Estimate_Channel3,Estimates3_grid] = Channel_Estimate_odd(Reference_Estimate3,numREperSubframe/2, CellRS_idx);
[Estimate_Channel4,Estimates4_grid] = Channel_Estimate_even(Reference_Estimate4,numREperSubframe/2, CellRS_idx2);

%% Recive Signal ZF Equalization & De-precoding

Equalize_PDSCH_A = zeros(length(PDSCH_idx),1); Equalize_PDSCH_B = zeros(length(PDSCH_idx),1);
Constant_A = zeros(length(PDSCH_idx),1);
for idx = 1:1:length(PDSCH_idx)
    CE_idx = PDSCH_idx(idx);
    Constant_A(idx,1) = Estimate_Channel1(CE_idx)*Estimate_Channel4(CE_idx) - Estimate_Channel2(CE_idx)*Estimate_Channel3(CE_idx);
    
    Equalize_PDSCH_A(idx,1) = (Estimate_Channel4(CE_idx)*Recive_LTE_Ant1(CE_idx) - Estimate_Channel2(CE_idx)*Recive_LTE_Ant2(CE_idx))/Constant_A(idx,1);
    Equalize_PDSCH_B(idx,1) = (Estimate_Channel1(CE_idx)*Recive_LTE_Ant2(CE_idx) - Estimate_Channel3(CE_idx)*Recive_LTE_Ant1(CE_idx))/Constant_A(idx,1);
end

% De-Precoding
De_Precoding = zeros(length(PDSCH_idx),2);

if (eNB.PDSCH.TxScheme == "CDD")
    for idx = 1:1:length(PDSCH_idx)/2
        De_Precoding(2*idx-1,1) = Equalize_PDSCH_A(2*idx-1) + Equalize_PDSCH_B(2*idx-1);
        De_Precoding(2*idx,1) = Equalize_PDSCH_A(2*idx) - Equalize_PDSCH_B(2*idx);
        
        De_Precoding(2*idx-1,2) = Equalize_PDSCH_A(2*idx-1) - Equalize_PDSCH_B(2*idx-1);
        De_Precoding(2*idx,2) = Equalize_PDSCH_A(2*idx) + Equalize_PDSCH_B(2*idx);
    end    
    
elseif (eNB.PDSCH.TxScheme == "SpatialMux")
    for idx = 1:1:length(PDSCH_idx)
        De_Precoding(idx,1) = Equalize_PDSCH_A(idx) + Equalize_PDSCH_B(idx);
        De_Precoding(idx,2) = Equalize_PDSCH_A(idx) - Equalize_PDSCH_B(idx);
    end    
end

Error_Deprecoding = sum(round(PDSCH_Data,4) ~= round(De_Precoding,4));

%% SSS ZF Equalization & De-precoding
if (Num_Subframe == 0)
    Equalize_SSS_A = zeros(length(PDSCH_idx),1); Equalize_SSS_B = zeros(length(PDSCH_idx),1);
    Constant_SSS  = zeros(length(SSS_idx)/2,1);
    for idx = 1:1:length(SSS_idx)
        CE_idx = SSS_idx(idx);
        Constant_SSS(idx,1) = Estimate_Channel1(CE_idx)*Estimate_Channel4(CE_idx) - Estimate_Channel2(CE_idx)*Estimate_Channel3(CE_idx);
        
        Equalize_SSS_A(idx,1) = (Estimate_Channel4(CE_idx)*Recive_LTE_Ant1(CE_idx) - Estimate_Channel2(CE_idx)*Recive_LTE_Ant2(CE_idx))/Constant_SSS(idx,1);
        Equalize_SSS_B(idx,1) = (Estimate_Channel1(CE_idx)*Recive_LTE_Ant2(CE_idx) - Estimate_Channel3(CE_idx)*Recive_LTE_Ant1(CE_idx))/Constant_SSS(idx,1);
    end    
    
    % De-Precoding
    De_Precoding_SSS = zeros(length(SSS_idx),2);
    if (eNB.PDSCH.TxScheme == "CDD")        
        for idx = 1:1:length(SSS_idx)/2
            De_Precoding_SSS(2*idx-1,1) = Equalize_SSS_A(2*idx-1) + Equalize_SSS_B(2*idx-1);
            De_Precoding_SSS(2*idx,1) = Equalize_SSS_A(2*idx) - Equalize_SSS_B(2*idx);
            
            De_Precoding_SSS(2*idx-1,2) = Equalize_SSS_A(2*idx-1) - Equalize_SSS_B(2*idx-1);
            De_Precoding_SSS(2*idx,2) = Equalize_SSS_A(2*idx) + Equalize_SSS_B(2*idx);
        end
    elseif (eNB.PDSCH.TxScheme == "SpatialMux")
        for idx = 1:1:length(SSS_idx)
            De_Precoding_SSS(idx,1) = Equalize_SSS_A(idx) + Equalize_SSS_B(idx);
            De_Precoding_SSS(idx,2) = Equalize_SSS_A(idx) - Equalize_SSS_B(idx);
        end
    end
end

%% LTE Transmition Mode 3(MIMO, 2 by 2, 10MHz, 64QAM) - PDSCH Decoding

PDSCH_Demodulation(:,1) = lteSymbolDemodulate(De_Precoding(:,1),MCS,'Soft');
PDSCH_Demodulation(:,2) = lteSymbolDemodulate(De_Precoding(:,2),MCS,'Soft');

PDSCH_Descramble = zeros(length(Descramble_Sequemce),2);
CRCDecode_Data = zeros(TransBlock_Size,2); Error_information = zeros(1,2);
for k = 1:1:2
    PDSCH_Descramble(:,k) = Descramble_Sequemce .* PDSCH_Demodulation(:,k);
    
    PDSCH_DerateMatch = lteRateRecoverTurbo(PDSCH_Descramble(:,k), TransBlock_Size, 0, eNB.PDSCH);
    
    % Turbo Decode
    Turbo_Decode_Data = lteTurboDecode(PDSCH_DerateMatch, 3);
    
    % Desegmentation
    [Desegment_Data, Desegment_err] = lteCodeBlockDesegment(Turbo_Decode_Data, length(CRC_insertion_Data));
    
    % CRC Decode
    [CRCDecode_Data(:,k), CRCDecode_err] = lteCRCDecode(Desegment_Data, '24A');
    
    Error_information(:,k) = sum(Information_bit(:,k) ~= CRCDecode_Data(:,k));
end