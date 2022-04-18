%% Initialization

clc
% clear;
% close all;

%% WIFI(802.11ac) Setting

CBW ='CBW20';
Scramble_Sequence = xlsread('IEEE802.11ac_Scramble_Sequence.xlsx');
Scramble_Sequence = Scramble_Sequence(93,:);

MCS = 8;
APEPLengthMCS = 1024;

if MCS == 0
    numBPSCS = 1; numCBPS = 52; numDBPS = 26;
    Rate = 0; Coding_Rate = 1/2; %Coding = '1/2';
elseif MCS == 1
    numBPSCS = 2; numCBPS = 104; numDBPS = 52; 
    Rate = 0; Coding_Rate = 1/2;  %Coding = '1/2';
elseif MCS == 2
    numBPSCS = 2; numCBPS = 104; numDBPS = 78;
    Rate = 2; Coding_Rate = 3/4;  %Coding = '3/4';
elseif MCS == 3
    numBPSCS = 4; numCBPS = 208; numDBPS = 104;
    Rate = 0; Coding_Rate = 1/2;  %Coding = '1/2';
elseif MCS == 4
    numBPSCS = 4; numCBPS = 208; numDBPS = 156; 
    Rate = 2; Coding_Rate = 3/4; %Coding = '3/4';
elseif MCS == 5
    numBPSCS = 6; numCBPS = 312; numDBPS = 208; 
    Rate = 1; Coding_Rate = 2/3; %Coding = '2/3';
elseif MCS == 6
    numBPSCS = 6; numCBPS = 312; numDBPS = 234; 
    Rate = 2; Coding_Rate = 3/4; %Coding = '3/4';
elseif MCS == 7
    numBPSCS = 6; numCBPS = 312; numDBPS = 260;
    Rate = 3; Coding_Rate = 5/6; %Coding = '5/6';
elseif MCS == 8
    numBPSCS = 8; numCBPS = 416; numDBPS = 312;
    Rate = 2; Coding_Rate = 3/4; %Coding = '3/4';
end

if (APEPLengthMCS ==0)
    Nsym = 0;
else
    Nsym = ceil((8*APEPLengthMCS+16+6)/numDBPS);
end
PHY_LENGTH = floor((Nsym*numDBPS-16-6)/8); % size = byte

%% Preamble Field

%% Preamble_VHT-SIG-B

VHT_SIG_B_LENGTH = ceil(APEPLengthMCS/4);
SIG_B_LENGTH = double(fliplr(dec2bin(VHT_SIG_B_LENGTH,17)))'-48;
SIG_B_LENGTH = [SIG_B_LENGTH; 1; 1; 1];
VHT_SIG_B_Tail = zeros(6,1);
VHT_SIG_B = [SIG_B_LENGTH; VHT_SIG_B_Tail];

% VHT_SIG_B_BCC
SIGB_L = length(VHT_SIG_B);
VHT_SIG_B_BCC = zeros(SIGB_L,2);
VHT_SIG_B_R = [VHT_SIG_B(end) VHT_SIG_B(end-1) VHT_SIG_B(end-2) VHT_SIG_B(end-3) VHT_SIG_B(end-4) VHT_SIG_B(end-5)];
for i = 1:SIGB_L
    VHT_SIG_B_BCC(i,1) = mod((VHT_SIG_B(i) + VHT_SIG_B_R(2) + VHT_SIG_B_R(3) + VHT_SIG_B_R(5) + VHT_SIG_B_R(6)), 2); 
    VHT_SIG_B_BCC(i,2) = mod((VHT_SIG_B(i) + VHT_SIG_B_R(1) + VHT_SIG_B_R(2) + VHT_SIG_B_R(3) + VHT_SIG_B_R(6)), 2); 
    VHT_SIG_B_R = circshift(VHT_SIG_B_R, [1,1]);
    VHT_SIG_B_R(1) = VHT_SIG_B(i);
end
VHT_SIG_B_BCCEncode = zeros(SIGB_L*2,1);
for n= 1:1:SIGB_L
    VHT_SIG_B_BCCEncode(2*n-1,1) = VHT_SIG_B_BCC(n,1);
    VHT_SIG_B_BCCEncode(2*n,1) = VHT_SIG_B_BCC(n,2);
end

numCBPSB = 52;
length_SIGB = (0:numCBPSB-1);
Interleaver_VHTSIGB = zeros(numCBPSB,1);

First_permutation = zeros(numCBPSB,1);
Second_permutation = zeros(numCBPSB,1);
for i = 1:1:1
    input_data = VHT_SIG_B_BCCEncode(1+numCBPSB*(i-1):numCBPSB*i);
    n = (numCBPSB/13)*mod(length_SIGB,13)+floor(length_SIGB/13);
    
    for m =1:1:length(length_SIGB)
        First_permutation(n(m)+1) = input_data(m);
    end
    
    n2 = floor(length_SIGB)+mod((length_SIGB+numCBPSB-floor(length_SIGB*13/numCBPSB)),1);
    
    for z = 1:1:length(length_SIGB)        
        Second_permutation(n2(z)+1) = First_permutation(z);  
    end
    Interleaver_VHTSIGB(:,i) = Second_permutation; 
end
Mapping_VHTSIGB = 2*Interleaver_VHTSIGB-1;
VHT_SIGNAL_B = [zeros(4,1); Mapping_VHTSIGB(1:7); 1; Mapping_VHTSIGB(8:20); 1; Mapping_VHTSIGB(21:26); 0; Mapping_VHTSIGB(27:32); 1; Mapping_VHTSIGB(33:45); -1; Mapping_VHTSIGB(46:52); zeros(3,1)];
Check_VHT_SIGNAL_B = [zeros(4,1); Mapping_VHTSIGB(1:7); 1; Mapping_VHTSIGB(8:20); 1; Mapping_VHTSIGB(21:26); 0; Mapping_VHTSIGB(27:32); 1; Mapping_VHTSIGB(33:45); -1; Mapping_VHTSIGB(46:52); zeros(3,1)];
% VHT_SIGNAL_B = ifft(circshift(VHT_SIGNAL_B',32)',64)*64/sqrt(56);
VHT_SIGNAL_B = ifft(circshift(VHT_SIGNAL_B',32)',64);
VHT_SIGNAL_B = [VHT_SIGNAL_B(49:64); VHT_SIGNAL_B];

%% Service field_Data field

ScramInit = zeros(7,1);
window_SIG_B = ones(8,1);
for n=1:length(SIG_B_LENGTH)
    a = xor(window_SIG_B(1),SIG_B_LENGTH(n));
    window_SIG_B(1:5) = window_SIG_B(2:6);
    window_SIG_B(6) = xor(window_SIG_B(7),a);
    window_SIG_B(7) = xor(window_SIG_B(end),a);
    window_SIG_B(end) = a;   
end
SIG_B_CRC = ~window_SIG_B;

Service_field = [ScramInit; 0; SIG_B_CRC];
Tail_bit = zeros(6,1);
Pad_Length = (Nsym*numDBPS)-(8*PHY_LENGTH)-22; 
Pad_bit = zeros(Pad_Length,1);

%% Data field

Data_information = randi([0 1],PHY_LENGTH*8,1);
Data_information_Test = zeros(PHY_LENGTH*8,1);
Test_Data = [1 0 0 1];
for sub = 1:1:2100
    for idx = 1:1:4
        Data_information_Test(idx+4*(sub-1)) = Test_Data(idx);
    end
end


% PSDU_plus = [Service_field; Data_information; Pad_bit; Tail_bit];
PSDU_plus = [Service_field; PSDU; Pad_bit; Tail_bit];

for sub = 1:1:42
    for idx = 1:1:200
%         DSP_Data_information(sub,idx) = Data_information(idx+200*(sub-1));
        DSP_Data_information(sub,idx) = PSDU(idx+200*(sub-1));
    end
end

%% Scramble_Data field

Scramble_data = zeros(length(PSDU_plus),1);
loop = floor(length(PSDU_plus)/127);
if (loop == 0)
        for n=1:1:(length(PSDU_plus))
            Scramble_data(n) = xor(PSDU_plus(n),Scramble_Sequence(n));
        end
else
    k = 0;
    while (k<loop)
        for n=1:1:127
            Scramble_data(n+127*k) = xor(PSDU_plus(n+127*k),Scramble_Sequence(n));
        end
        k = k+1;
    end
    
    for n = 1:1:(length(PSDU_plus)-127*loop)
        Scramble_data(n+127*loop) = xor(PSDU_plus(n+127*loop),Scramble_Sequence(n));
    end
end
% Error_Scramble =sum(Check_Scramble~=Scramble_data);
Scramble_data =[Scramble_data(1:(length(PSDU_plus)-6)); Tail_bit];

%% Binary Convolution Code_Data field

PSDU_L = length(PSDU_plus);
VHT_Data_BCC = zeros(PSDU_L,2);
Data_R = [Scramble_data(end) Scramble_data(end-1) Scramble_data(end-2) Scramble_data(end-3) Scramble_data(end-4) Scramble_data(end-5)];
for i = 1:PSDU_L
    VHT_Data_BCC(i,1) = mod((Scramble_data(i) + Data_R(2) + Data_R(3) + Data_R(5) + Data_R(6)), 2); 
    VHT_Data_BCC(i,2) = mod((Scramble_data(i) + Data_R(1) + Data_R(2) + Data_R(3) + Data_R(6)), 2); 
    Data_R = circshift(Data_R, [1,1]);
    Data_R(1) = Scramble_data(i);
end

%% Binary Convolution Code_Coding Rate

Coding_Rate_length = length(PSDU_plus)*(1/Coding_Rate);
VHT_Data_BCCEncode = zeros(Coding_Rate_length,1);

switch Rate
    case 0 % 1/2
        for n = 1:1:PSDU_L
            VHT_Data_BCCEncode(2*n-1,1) = VHT_Data_BCC(n,1);
            VHT_Data_BCCEncode(2*n,1) = VHT_Data_BCC(n,2);
        end
    case 1 % 2/3
        for n = 1:1:(Coding_Rate_length/3)
            VHT_Data_BCCEncode(3*n-2) = VHT_Data_BCC(2*n-1,1);
            VHT_Data_BCCEncode(3*n-1) = VHT_Data_BCC(2*n-1,2);
            VHT_Data_BCCEncode(3*n) = VHT_Data_BCC(2*n,1);
        end
    case 2 % 3/4
        for n = 1:1:(Coding_Rate_length/4)
            VHT_Data_BCCEncode(4*n-3) = VHT_Data_BCC(3*n-2,1);            
            VHT_Data_BCCEncode(4*n-2) = VHT_Data_BCC(3*n-2,2);            
            VHT_Data_BCCEncode(4*n-1) = VHT_Data_BCC(3*n-1,1);
            VHT_Data_BCCEncode(4*n) = VHT_Data_BCC(3*n,2);            
        end
    case 3 % 5/6
        for n = 1:1:(Coding_Rate_length/12)
            VHT_Data_BCCEncode(12*n-11) = VHT_Data_BCC(10*n-9,1);
            VHT_Data_BCCEncode(12*n-9) = VHT_Data_BCC(10*n-8,1);
            VHT_Data_BCCEncode(12*n-7) = VHT_Data_BCC(10*n-6,1);            
            VHT_Data_BCCEncode(12*n-5) = VHT_Data_BCC(10*n-4,1);
            VHT_Data_BCCEncode(12*n-3) = VHT_Data_BCC(10*n-3,1);
            VHT_Data_BCCEncode(12*n-1) = VHT_Data_BCC(10*n-1,1);
            
            VHT_Data_BCCEncode(12*n-10) = VHT_Data_BCC(10*n-9,2);
            VHT_Data_BCCEncode(12*n-8) = VHT_Data_BCC(10*n-7,2);
            VHT_Data_BCCEncode(12*n-6) = VHT_Data_BCC(10*n-5,2);
            VHT_Data_BCCEncode(12*n-4) = VHT_Data_BCC(10*n-4,2);
            VHT_Data_BCCEncode(12*n-2) = VHT_Data_BCC(10*n-2,2);
            VHT_Data_BCCEncode(12*n) = VHT_Data_BCC(10*n,2);            
        end
end

%% Binary Convolution Code Interleave_Data field

Interleave_length = (0:numCBPS-1);
Interleaver_data = zeros(numCBPS,Nsym);
s = max(1,numBPSCS/2);

First_permutation = zeros(numCBPS,1);
Second_permutation = zeros(numCBPS,1);

for i = 1:1:Nsym
    input_data = VHT_Data_BCCEncode(1+numCBPS*(i-1):numCBPS*i);
    n = (numCBPS/13)*mod(Interleave_length,13)+floor(Interleave_length/13);
    
    for m =1:1:length(Interleave_length)
        First_permutation(n(m)+1) = input_data(m);
    end
    n2 = s*floor(Interleave_length/s)+mod((Interleave_length+numCBPS-floor(Interleave_length*13/numCBPS)),s);
 
    for z = 1:1:length(Interleave_length)
        Second_permutation(n2(z)+1) = First_permutation(z);
    end
    Interleaver_data(:,i) = Second_permutation;    
end

Data = zeros(length(VHT_Data_BCCEncode),1);
for k = 0:1: Nsym-1
    for m = 1:1:numCBPS
        Data(m+k*numCBPS) = Interleaver_data(m,k+1);
    end
end

%% Modulation

Modulation_length = Coding_Rate_length/numBPSCS;
Modulation = zeros(Modulation_length,2);
Input_DATA = 2*Data-1;

for i=1:1:Modulation_length
    m1 = Input_DATA(8*i-4)-2;
    m1 = abs(m1)*Input_DATA(8*i-5)-4;
    m1 = abs(m1)*Input_DATA(8*i-6)-8;
    Modulation(i,1) = abs(m1)*Input_DATA(8*i-7);
    
    m2 = Input_DATA(8*i)-2;
    m2 = abs(m2)*Input_DATA(8*i-1)-4;
    m2 = abs(m2)*Input_DATA(8*i-2)-8;
    Modulation(i,2) = abs(m2)*Input_DATA(8*i-3);
end

% Modulation_Data = (Modulation(:,1) + 1j*Modulation(:,2));
Modulation_Data = (Modulation(:,1) + 1j*Modulation(:,2))/sqrt(170);
% Modulation_Data = (Modulation(:,1) + 1j*Modulation(:,2))/13;
Modulation_Data_real = real(Modulation_Data * 2^12);
Modulation_Data_imag = imag(Modulation_Data * 2^12);

%% Zero_Pading, Pilot Insertion_Data field
Polarity = [1;1;1;1; -1;-1;-1;1; -1;-1;-1;-1; 1;1;-1;1; -1;-1;1;1; -1;1;1;-1; 1;1;1;1; 1;1;-1;1;
1;1;-1;1; 1;-1;-1;1; 1;1;-1;1; -1;-1;-1;1; -1;1;-1;-1; 1;-1;-1;1; 1;1;1;1; -1;-1;1;1;
-1;-1;1;-1; 1;-1;1;1; -1;-1;-1;1; 1;-1;-1;-1; -1;1;-1;-1; 1;-1;1;1; 1;1;-1;1; -1;1;-1;1;
-1;-1;-1;-1; -1;1;-1;1; 1;-1;1;-1; 1;1;1;-1; -1;1;-1;-1; -1;1;1;1; -1;-1;-1;-1; -1;-1;-1]; % Num = 127

Pilot = [1, 1, 1, -1; 1, 1, -1, 1; 1, -1, 1, 1; -1, 1, 1, 1];
Pilot = Pilot';
length_Transmit_Data = length(Modulation_Data)+(12*Nsym);
Transmit_Data = zeros(length_Transmit_Data,1);

for n = 1:1:Nsym
    for i = 1:1:7
        Transmit_Data(i+4+64*(n-1)) = Modulation_Data(i+52*(n-1));
        Transmit_Data(i+54+64*(n-1)) = Modulation_Data(i+45+52*(n-1));
    end
    for i = 1:1:13
        Transmit_Data(i+12+64*(n-1)) = Modulation_Data(i+7+52*(n-1));
        Transmit_Data(i+40+64*(n-1)) = Modulation_Data(i+32+52*(n-1));
    end
    for i = 1:1:6
        Transmit_Data(i+26+64*(n-1)) = Modulation_Data(i+20+52*(n-1));
        Transmit_Data(i+33+64*(n-1)) = Modulation_Data(i+26+52*(n-1));
    end
    
    if n<124
        Polarity_Pilot = Pilot(:,mod(n-1,4)+1)*Polarity(n+4);
    else
        Polarity_Pilot = Pilot(:,mod(n-1,4)+1)*Polarity((mod(n-1,123)+1));
    end
    
    for i = 1:1:4
        Transmit_Data(14*i-2+64*(n-1)) = Polarity_Pilot(i);
    end
end

IFFT_input = zeros(64,Nsym); IFFT_ = zeros(64,Nsym); IFFT_Data = zeros(80,Nsym); 
% for k = 1:1:Nsym
%     IFFT_input(:,k) = Transmit_Data(1+64*(k-1):1:64+64*(k-1));
%     IFFT_(:,k) = ifft(IFFT_input(:,k),64)*64 / sqrt(56);
%     IFFT_Data(:,k) = [IFFT_(49:64,k); IFFT_(:,k)];
% end
Transmit_Data_row = reshape(Transmit_Data,64,27);
vhtdata_ifft =ifft(ifftshift(Transmit_Data_row,1));

% vhtdata_ifft =(1/sqrt(56))* 64 * ifft(ifftshift(Transmit_Data_row,1));



vhtdata_ifft_cp = zeros(80,Nsym);
for cpnn = 1:Nsym
    vhtdata_ifft_cp(:,cpnn) = [vhtdata_ifft(49:64,cpnn);vhtdata_ifft(:,cpnn)];
end
VHT_Data = reshape(vhtdata_ifft_cp,[80*Nsym,1]);

% VHT_Data = zeros(80*Nsym,1);
% for k=1:1:Nsym
%    for i = 1:1:80
%        VHT_Data(i+80*(k-1),1) = IFFT_Data(i,k);
%    end    
% end

real_part = real(VHT_Data*2^9);
imag_part = imag(VHT_Data*2^9);

real_Wave = fi(int16(real(VHT_Data*2^12)));
real_Wave_max = max(real_Wave);
real_Wave_min = min(real_Wave);

%% Preable Field Encoding Processing
%% Preamble Legacy Signal field

L_SIG_Rate = [1; 1; 0; 1];
L_SIG_Tail = zeros(6,1);

T_Data = 4*(Nsym);
TXTIME = 40 + T_Data;
L_LENGTH = ceil((TXTIME-20)/4)*3-3;
L_LENGTH = double(fliplr(dec2bin(L_LENGTH,12)))'-48;
L_SIG = [L_SIG_Rate; 0; L_LENGTH; L_SIG_Tail];

if (rem(sum(L_SIG),2)==0)
    L_SIG_Parity = 0;
else
    L_SIG_Parity = 1;
end
L_SIG = [L_SIG_Rate; 0; L_LENGTH; L_SIG_Parity; L_SIG_Tail];

%% Legacy_SIG BCC
SIG_L = length(L_SIG);
Legacy_SIG_BCC = zeros(SIG_L,2);
LSIG_R = [L_SIG(end) L_SIG(end-1) L_SIG(end-2) L_SIG(end-3) L_SIG(end-4) L_SIG(end-5)];
for i = 1:1:SIG_L
    Legacy_SIG_BCC(i,1) = mod((L_SIG(i) + LSIG_R(2) + LSIG_R(3) + LSIG_R(5) + LSIG_R(6)), 2); 
    Legacy_SIG_BCC(i,2) = mod((L_SIG(i) + LSIG_R(1) + LSIG_R(2) + LSIG_R(3) + LSIG_R(6)), 2); 
    LSIG_R = circshift(LSIG_R, [1,1]);
    LSIG_R(1) = L_SIG(i);
end

Legacy_SIG_BCCEncode = zeros(2*SIG_L,1);
for n= 1:1:SIG_L
    Legacy_SIG_BCCEncode(2*n-1,1) = Legacy_SIG_BCC(n,1);
    Legacy_SIG_BCCEncode(2*n,1) = Legacy_SIG_BCC(n,2);
end

numCBPSS = 48;
length_LSIG = (0:numCBPSS-1);
Interleaver_LSIG = zeros(numCBPSS,1);
s = 1;

First_permutation = zeros(numCBPSS,1);
Second_permutation = zeros(numCBPSS,1);

for i = 1:1:1
    input_data = Legacy_SIG_BCCEncode(1+numCBPSS*(i-1):numCBPSS*i)';
    n = (numCBPSS/16)*mod(length_LSIG,16)+floor(length_LSIG/16);
    
    for m =1:1:length(length_LSIG)
        First_permutation(n(m)+1) = input_data(m);
    end
    n2 = s*floor(length_LSIG/s)+mod((length_LSIG+numCBPSS-floor(length_LSIG*16/numCBPSS)),s);
 
    for z = 1:1:length(length_LSIG)
        Second_permutation(n2(z)+1) = First_permutation(z);
    end
    Interleaver_LSIG(:,i) = Second_permutation;       
end
Mapping_LSIG = 2*Interleaver_LSIG-1;
Legacy_SIG = [zeros(6,1); Mapping_LSIG(1:5);1;Mapping_LSIG(6:18);1;Mapping_LSIG(19:24);0;Mapping_LSIG(25:30);1;Mapping_LSIG(31:43);-1;Mapping_LSIG(44:48); zeros(5,1)];
% Legacy_SIG_ifft = ifft(circshift(Legacy_SIG,32),64)*64/sqrt(52);
Legacy_SIG_ifft = ifft(circshift(Legacy_SIG,32),64);
Legacy_SIGNAL = [Legacy_SIG_ifft(49:64); Legacy_SIG_ifft];

%% Preamble field (VHT Field)
%% VHT Signal A
BW = 20;
if BW ==20
    BW = [0; 0];
elseif BW == 40
    BW = [0; 1];
elseif BW == 80
    BW = [1; 0];
elseif BW == 160
    BW = [1; 1];
end
GroupID = 63;
STBC = 0; NSTS = 0; TXOP = 0;
PartialAID = 275;

GroupID = double(fliplr(dec2bin(GroupID,6)))'-48;
NSTS = double(fliplr(dec2bin(NSTS,3)))'-48;
PartialAID = double(fliplr(dec2bin(PartialAID,9)))'-48;

VHT_SIG_A1 = [BW; 1; STBC; GroupID; NSTS; PartialAID; TXOP; 1];

ShortGI = [0; 0];
Coding = 0; % 0_BCC, 1_LDPC
LDPC = 0; Beamform = 0;
MCS = double(fliplr(dec2bin(MCS,4)))'-48;
VHT_SIG_A_Tail = zeros(6,1);
VHT_SIG_A2 = [ShortGI; Coding; LDPC; MCS; Beamform; 1];
VHT_SIG_A = [VHT_SIG_A1; VHT_SIG_A2];

% VHT_SIG_A_CRC
window_SIG_A = ones(8,1);
for n=1:length(VHT_SIG_A)
    a = xor(window_SIG_A(1),VHT_SIG_A(n));
    window_SIG_A(1:5) = window_SIG_A(2:6);
    window_SIG_A(6) = xor(window_SIG_A(7),a);
    window_SIG_A(7) = xor(window_SIG_A(end),a);
    window_SIG_A(end) = a;   
end
VHT_SIG_A_CRC = ~window_SIG_A;

VHT_SIG_A = [VHT_SIG_A1; VHT_SIG_A2; VHT_SIG_A_CRC; VHT_SIG_A_Tail];

%% VHT-SIG-A BCC
SIGA_L = length(VHT_SIG_A);
VHT_SIG_A_BCC = zeros(SIGA_L,2);
SIGA_R = [VHT_SIG_A(end) VHT_SIG_A(end-1) VHT_SIG_A(end-2) VHT_SIG_A(end-3) VHT_SIG_A(end-4) VHT_SIG_A(end-5)];
for i = 1:SIGA_L
    VHT_SIG_A_BCC(i,1) = mod((VHT_SIG_A(i) + SIGA_R(2) + SIGA_R(3) + SIGA_R(5) + SIGA_R(6)), 2); 
    VHT_SIG_A_BCC(i,2) = mod((VHT_SIG_A(i) + SIGA_R(1) + SIGA_R(2) + SIGA_R(3) + SIGA_R(6)), 2); 
    SIGA_R = circshift(SIGA_R, [1,1]);
    SIGA_R(1) = VHT_SIG_A(i);
end
VHT_SIG_A_BCCEncode = zeros(SIGA_L*2,1);
for n= 1:1:SIGA_L
    VHT_SIG_A_BCCEncode(2*n-1,1) = VHT_SIG_A_BCC(n,1);
    VHT_SIG_A_BCCEncode(2*n,1) = VHT_SIG_A_BCC(n,2);
end

length_SIGA = (0:numCBPSS-1);
Interleaver_VHTSIGA = zeros(numCBPSS,1);

First_permutation = zeros(numCBPSS,1);
Second_permutation = zeros(numCBPSS,1);
for i = 1:1:2
    input_data = VHT_SIG_A_BCCEncode(1+numCBPSS*(i-1):numCBPSS*i);   
    n = (numCBPSS/16)*mod(length_SIGA,16)+floor(length_SIGA/16);
    
    for m =1:1:length(length_SIGA)
        First_permutation(n(m)+1) = input_data(m);
    end
    n2 = s*floor(length_SIGA/s)+mod((length_SIGA+numCBPSS-floor(length_SIGA*16/numCBPSS)),s);
 
    for z = 1:1:length(length_SIGA)
        Second_permutation(n2(z)+1) = First_permutation(z);
    end
    Interleaver_VHTSIGA(:,i) = Second_permutation;    
end

Mapping_VHTSIGA(:,1) = 2*Interleaver_VHTSIGA(:,1)-1;
Mapping_VHTSIGA(:,2) = (2*Interleaver_VHTSIGA(:,2)-1);

VHT_SIG_A1 = [zeros(6,1); Mapping_VHTSIGA(1:5,1);1;Mapping_VHTSIGA(6:18,1);1;Mapping_VHTSIGA(19:24,1);0;Mapping_VHTSIGA(25:30,1);1;Mapping_VHTSIGA(31:43,1);-1;Mapping_VHTSIGA(44:48,1); zeros(5,1)];
VHTA1= [zeros(6,1); Mapping_VHTSIGA(1:5,1);1;Mapping_VHTSIGA(6:18,1);1;Mapping_VHTSIGA(19:24,1);0;Mapping_VHTSIGA(25:30,1);1;Mapping_VHTSIGA(31:43,1);-1;Mapping_VHTSIGA(44:48,1); zeros(5,1)];
VHT_SIG_A2 = [zeros(6,1); Mapping_VHTSIGA(1:5,2);1;Mapping_VHTSIGA(6:18,2);1;Mapping_VHTSIGA(19:24,2);0;Mapping_VHTSIGA(25:30,2);1;Mapping_VHTSIGA(31:43,2);-1;Mapping_VHTSIGA(44:48,2); zeros(5,1)];
VHTA2 =[zeros(6,1); Mapping_VHTSIGA(1:5,2);1;Mapping_VHTSIGA(6:18,2);1;Mapping_VHTSIGA(19:24,2);0;Mapping_VHTSIGA(25:30,2);1;Mapping_VHTSIGA(31:43,2);-1;Mapping_VHTSIGA(44:48,2); zeros(5,1)];
VHT_SIG_A2 = VHT_SIG_A2*exp(1j*pi/2);

% VHT_SIG_A1 = ifft(circshift(VHT_SIG_A1',32)',64)*(64/sqrt(52));
% VHT_SIG_A2 = ifft(circshift(VHT_SIG_A2',32)',64)*(64/sqrt(52));
VHT_SIG_A1 = ifft(circshift(VHT_SIG_A1',32)',64);
VHT_SIG_A2 = ifft(circshift(VHT_SIG_A2',32)',64);
Preamble_VHT_SIG_A = [VHT_SIG_A1(49:64); VHT_SIG_A1; VHT_SIG_A2(49:64); VHT_SIG_A2];

%%
L_STF_seq = [0,0,0,0,0,0,0,0,1+1j,0,0,0,-1-1j,0,0,0,1+1j,0,0,0,-1-1j,0,0,0,-1-1j,0,0,0,1+1j,0,0,0,0,0,0,0,-1-1j,0,0,0,-1-1j,0,0,0,1+1j,0,0,0,1+1j,0,0,0,1+1j,0,0,0,1+1j,0,0,0,0,0,0,0];

L_STF_seq = circshift(L_STF_seq,32);
% L_STF_seq = circshift(L_STF_seq,32)*sqrt(1/2);
% L_STF(:,1) = ifft(L_STF_seq,64)*64/sqrt(12);
L_STF(:,1) = ifft(L_STF_seq,64);
Legacy_STF = [L_STF(49:64); L_STF; L_STF(49:64); L_STF];
Preamble_VHT_STF = [L_STF(49:64); L_STF];
Preamble_VHT_STF_scale = [L_STF(49:64); L_STF];


L_LTF_seq = [0,0,0,0,0,0,1,1,-1,-1,1,1,-1,1,-1,1,1,1,1,1,1,-1,-1,1,1,-1,1,-1,1,1,1,1,0,1,-1,-1,1,1,-1,1,-1,1,-1,-1,-1,-1,-1,1,1,-1,-1,1,-1,1,-1,1,1,1,1,0,0,0,0,0];
L_LTF_seq = circshift(L_LTF_seq,32);
% L_LTF(:,1) = ifft(L_LTF_seq,64)*(64/sqrt(52));
L_LTF(:,1) = ifft(L_LTF_seq,64);
Legacy_LTF = [L_LTF(33:64); L_LTF; L_LTF];
% Legacy_LTF = [L_LTF(49:64); L_LTF; L_LTF(49:64); L_LTF];
Legacy_LTF_scale = Legacy_LTF;

VHT_LTF_seq = [0 0 0 0 1 1 1 1 -1 -1 1 1 -1 1 -1 1 1 1 1 1 1 -1 -1 1 1 -1 1 -1 1 1 1 1 0 1 -1 -1 1 1 -1 1 -1 1 -1 -1 -1 -1 -1 1 1 -1 -1 1 -1 1 -1 1 1 1 1 -1 -1 0 0 0];
VHT_LTF_seq = circshift(VHT_LTF_seq,32);
% VHT_LTF(:,1) = ifft(VHT_LTF_seq,64)*(64/sqrt(56));
VHT_LTF(:,1) = ifft(VHT_LTF_seq,64);
Preamble_VHT_LTF = [VHT_LTF(49:64); VHT_LTF];
Preamble_VHT_LTF_scale = Preamble_VHT_LTF;

Preamble_Waveform = [Legacy_STF; Legacy_LTF_scale; Legacy_SIGNAL; Preamble_VHT_SIG_A; Preamble_VHT_STF; Preamble_VHT_LTF_scale; VHT_SIGNAL_B];
Preamble_length = length(Preamble_Waveform)/80;

Preamble_Waveform_scale = Preamble_Waveform * 2^16;
VHT_Data_scale = VHT_Data * 2^16;
Antenna_Power = 2^0;
% WIFI_Waveform = [Preamble_Waveform; VHT_Data] * Antenna_Power;
WIFI_Waveform = [Preamble_Waveform_scale; VHT_Data_scale] / Antenna_Power;
Check_Waveform = WIFI_Waveform/Antenna_Power;
figure(3)
plot(abs(WIFI_Waveform))


figure(18)
plot(abs(WIFI_Waveform))
hold on
plot(abs(RX_Data))


Scaling = 20;
Real_Part = real(WIFI_Waveform);
Imag_Part = imag(WIFI_Waveform);
frameL = [Imag_Part,Real_Part];
frameLL = reshape(frameL',5920,1);

[startOffset,meric] = wlanSymbolTimingEstimate(WIFI_Waveform,'CBW20'); %% rxSig  RX_Data WIFI_Waveform
% RX_Data = RX_Data(startOffset+1:startOffset+2960);
figure(22)
plot(meric)

%% Fading Channel
% sampleRate20MHz = 20e6;
% maxDopplerShift  = 60;
% delayVector = (0:5:15)*1e-6; % Discrete delays of four-path channel (s)
% % gainVector  = [0 -3 -6 -9];  % Average path gains (dB)
% gainVector  = [0 -2 -4 -6];
% KFactor = 10;            % Linear ratio of specular power to diffuse power
% specDopplerShift = 100;  % Doppler shift of specular component (Hz)
% 
% rayChan = comm.RayleighChannel( ...
%     'SampleRate',          sampleRate20MHz, ...
%     'PathDelays',          delayVector, ...
%     'AveragePathGains',    gainVector, ...
%     'MaximumDopplerShift', maxDopplerShift, ...
%     'RandomStream',        'mt19937ar with seed', ...
%     'Seed',                10, ...
%     'PathGainsOutputPort', true);
% 
% EbN0_dB = 0:2:16;
% EbN0=10.^(EbN0_dB(1)/10);
% sigma = sqrt((1/2)*(1/(EbN0))*(1/log2(256)));
% Noise = (randn(1,length(WIFI_Waveform))+1j*randn(1,length(WIFI_Waveform)));
% AWGN(:,1) = sigma*Noise/sqrt(2);
% 
% WIFI_Waveform = rayChan(WIFI_Waveform) + AWGN;
% Tx_Waveform = WIFI_Waveform;
% 
% % WIFI_Waveform = [WIFI_Waveform(801:2960); WIFI_Waveform(801:2960); WIFI_Waveform];
% 
% %% Sync - L-STF
% 
% Stand_STF = L_STF;
% STF_corr = zeros(length(WIFI_Waveform),1);
% for i=1:1:length(WIFI_Waveform)
%     if (i<length(Stand_STF))
%         STF_corr(i,1) = sum(WIFI_Waveform(1:i).*conj(Stand_STF(length(Stand_STF)+1-i:length(Stand_STF))));
%     else
%         STF_corr(i,1) = sum(WIFI_Waveform(i-(length(Stand_STF)-1):i).*conj(Stand_STF));        
%     end
% end
% 
% grid_STF = 1;
% while 1
%     Value = abs(STF_corr(grid_STF))-50*Antenna_Power+2;
%     if (Value >0)
%         break
%     end
%     grid_STF = grid_STF + 1;
% end
% 
% Sync_Data = WIFI_Waveform(grid_STF-64+1:end);
% 
% % figure(1)
% % plot(abs(STF_corr)-50*Antenna_Power+2)
% % grid
% % xlabel('Time Domain Sample Number')
% % ylabel('Magnitude [dB]')
% % axis([0 2980 0 30]);
% 
% %% Sync - L-LTF
% 
% Stand_LTF = L_LTF;
% LTF_corr = zeros(length(Sync_Data),1);
% for i=1:1:length(Sync_Data)
%     if (i<length(Stand_LTF))
%         LTF_corr(i,1) = sum(Sync_Data(1:i).*conj(Stand_LTF(length(Stand_LTF)+1-i:length(Stand_LTF))));
%     else
%         LTF_corr(i,1) = sum(Sync_Data(i-(length(Stand_LTF)-1):i).*conj(Stand_LTF));        
%     end
% end
% 
% grid_LTF = 1;
% while (grid_LTF<length(LTF_corr))
%     Value = abs(LTF_corr(grid_LTF))-45*Antenna_Power+2;
%     if (Value >0)
%         break
%     end
%     grid_LTF = grid_LTF + 1;
% end
% 
% % figure(2)
% % plot(abs(LTF_corr)-45*Antenna_Power+2)
% % grid
% % xlabel('Time Domain Sample Number')
% % ylabel('Magnitude [dB]')
% % axis([0 2980 0 10]);
% 
% %% Final Sync
% 
% if (grid_LTF == 320 || grid_LTF == 256)
%     Waveform = Sync_Data;
%     
% elseif (grid_LTF == 304)
%         Sync_Data = WIFI_Waveform(grid_STF-64+1-16:end);
%         Waveform = Sync_Data;    
%     
% elseif (grid_LTF == 288)
%         Sync_Data = WIFI_Waveform(grid_STF-64+1-32:end);
%         Waveform = Sync_Data;
%     
% elseif (grid_LTF == 272)
%         Sync_Data = WIFI_Waveform(grid_STF-64+1-48:end);
%         Waveform = Sync_Data;   
% 
% elseif (grid_LTF == length(LTF_corr))
%     Sync_Data = WIFI_Waveform(grid_STF-624+1:end);
%     
%     Stand_LTF = L_LTF;
%     LTF_corr = zeros(length(Sync_Data),1);
%     for i=1:1:length(Sync_Data)
%         if (i<length(Stand_LTF))
%             LTF_corr(i,1) = sum(Sync_Data(1:i).*conj(Stand_LTF(length(Stand_LTF)+1-i:length(Stand_LTF))));
%         else
%             LTF_corr(i,1) = sum(Sync_Data(i-(length(Stand_LTF)-1):i).*conj(Stand_LTF));
%         end
%     end
%     
%     grid_LTF = 1;
%     while (grid_LTF<length(LTF_corr))
%         Value = abs(LTF_corr(grid_LTF))-45*Antenna_Power+2;
%         if (Value >0)
%             break
%         end
%         grid_LTF = grid_LTF + 1;
%     end
%     
%     if (grid_LTF == 320 || grid_LTF == 256)
%         Waveform = Sync_Data;
%     elseif (grid_LTF == 304)
%         Sync_Data = WIFI_Waveform(grid_STF-624+1-16:end);
%         Waveform = Sync_Data;
%     elseif (grid_LTF == 288)
%         Sync_Data = WIFI_Waveform(grid_STF-624+1-32:end);
%         Waveform = Sync_Data;
%     elseif (grid_LTF == 272)
%         Sync_Data = WIFI_Waveform(grid_STF-624+1-48:end);
%         Waveform = Sync_Data;
%     end
% end
% 
% % figure(3)
% % plot(abs(LTF_corr)-45*Antenna_Power+2)
% % grid
% % xlabel('Time Domain Sample Number')
% % ylabel('Magnitude [dB]')
% % axis([0 2980 0 3]);
% 
% %% Decoding Processing
% 
% % CP Remove, FFT, Shift
% Symbol = Preamble_length + Nsym;
% FFT_Input = zeros(64,Symbol); FFT_ = zeros(64,Symbol);
% Detransmit_ = zeros(64,Symbol); Detransmit_Data = zeros(64*Symbol,1);
% for k = 1:1:Symbol
%         FFT_Input(:,k) = Waveform(17+80*(k-1):1:80+80*(k-1));
% %         FFT_Input(:,k) = rxwaveform(17+80*(k-1):1:80+80*(k-1));
%         FFT_(:,k) = fft(FFT_Input(:,k),64);
%         Detransmit_(:,k) = circshift(FFT_(:,k),32);
%         for i = 1:1:64
%             Detransmit_Data(i+64*(k-1)) = Detransmit_(i,k);
%         end
% end
% Detransmit_Data = round(Detransmit_Data,3);
% 
% %% Channel Estimation & Compensation
% 
% Pilot_Research = ones(64*Symbol,1);
% Detransmit_Pilot = zeros(64*Symbol,1);
% for n = 1:1:Nsym    
%     if n<124
%         Polarity_Pilot = Pilot(:,mod(n-1,4)+1)*Polarity(n+4);
%     else
%         Polarity_Pilot = Pilot(:,mod(n-1,4)+1)*Polarity((mod(n-1,123)+1));
%     end
%     
%     for i = 1:1:4
%         Detransmit_Pilot(14*i-2+64*(n-1)+640,1) = Detransmit_Data(14*i-2+64*(n-1)+640,1);
%         Pilot_Research(14*i-2+64*(n-1)+640,1) = Polarity_Pilot(i);
%     end
% end
% 
% Pattern = [1 1 1 -1];
% for n = 1:1:3
%     for i = 1:1:4
%          Detransmit_Pilot(14*i-2+64*(n-1)+256,1) = Detransmit_Data(14*i-2+64*(n-1)+256,1);
%          Pilot_Research(14*i-2+64*(n-1)+256,1) = Pattern(i);
%          Detransmit_Pilot(14*i-2+64*(6-1)+256,1) = Detransmit_Data(14*i-2+64*(6-1)+256,1);
%          Pilot_Research(14*i-2+64*(6-1)+256,1) = Pattern(i);
%     end    
% end
% Estimate_Channel =  (Detransmit_Pilot ./ Pilot_Research);
% 
% %% Channel Estimate : Interpolation & Extrapolation
% 
% for i = 1:1:256
%     Estimate_Channel(i) = 1;
% end
% for i = 1:1:128
%     Estimate_Channel(i+448) = 1;
% end
% 
% % Leagacy field Signal, VHT Signal A1, A2
% for n = 1:1:3
%     Inclination_A = (Estimate_Channel(282+64*(n-1),1)-Estimate_Channel(268+64*(n-1),1))/14;
%     for i = 1:1:26
%         Estimate_Channel(i+256+64*(n-1),1) = Inclination_A*i - Inclination_A*1 +Estimate_Channel(268+64*(n-1),1);
%     end    
%     Inclination_B = (Estimate_Channel(296+64*(n-1),1)-Estimate_Channel(282+64*(n-1),1))/14;
%     for i = 1:1:15
%         Estimate_Channel(i+281+64*(n-1),1) = Inclination_B*i - Inclination_B*1 +Estimate_Channel(282+64*(n-1),1);
%     end    
%     Inclination_C = (Estimate_Channel(310+64*(n-1),1)-Estimate_Channel(296+64*(n-1),1))/14;
%     for i = 1:1:25
%         Estimate_Channel(i+295+64*(n-1),1) = Inclination_C*i - Inclination_C*1 +Estimate_Channel(296+64*(n-1),1);
%     end
% end
% 
% % VHT Signal B
% Inclination_A = (Estimate_Channel(602,1)-Estimate_Channel(588,1))/14;
% for i = 1:1:26
%     Estimate_Channel(i+576,1) = Inclination_A*i - Inclination_A*1 +Estimate_Channel(588,1);
% end
% Inclination_B = (Estimate_Channel(616,1)-Estimate_Channel(602,1))/14;
% for i = 1:1:15
%     Estimate_Channel(i+601,1) = Inclination_B*i - Inclination_B*1 +Estimate_Channel(602,1);
% end
% Inclination_C = (Estimate_Channel(630,1)-Estimate_Channel(616,1))/14;
% for i = 1:1:25
%     Estimate_Channel(i+615,1) = Inclination_C*i - Inclination_C*1 +Estimate_Channel(616,1);
% end
% 
% % % Date field
% % for n = 1:1: Nsym
% % 
% % end
% 
% Estimate_Data = Detransmit_Data./Estimate_Channel;
% 
% %% Preamble Field Decoding
% 
% Deinterleave_Table = xlsread('Deinterleave_Table.xlsx');
% trellis = poly2trellis(7, [133 171]);
% 
% %% Leagacy Filed Signal Decoding Processing
% 
% DeLeagacy_SIG = Estimate_Data(257:320);
% for i = 1:1:64
%    if (DeLeagacy_SIG(i) <=0)
%        DeLeagacy_SIG(i) = 0;
%    else
%        DeLeagacy_SIG(i) = 1;
%    end
% end
% 
% DeLeagacy = zeros(48,1);
% for i = 1:1:5
%     DeLeagacy(i) = DeLeagacy_SIG(i+6);
%     DeLeagacy(i+43) = DeLeagacy_SIG(i+54);
% end
% for i = 1:1:13
%     DeLeagacy(i+5) = DeLeagacy_SIG(i+12);
%     DeLeagacy(i+30) = DeLeagacy_SIG(i+40);
% end
% for i = 1:1:6
%     DeLeagacy(i+24) = DeLeagacy_SIG(i+33);
%     DeLeagacy(i+18) = DeLeagacy_SIG(i+26);
% end
% Error_LSIG = sum(Interleaver_LSIG ~= DeLeagacy);
% LSIG_DeinterleaveTable = Deinterleave_Table(1:48,2)+1;
% 
% LSIG_Deinterleave = zeros(48,1);
% for i = 1:1:48
%     n = LSIG_DeinterleaveTable(i);
%     LSIG_Deinterleave(n,1) = DeLeagacy(i);
% end
% 
% Error_LSIG_Deinterleave = sum(LSIG_Deinterleave~=Legacy_SIG_BCCEncode);
% 
% Viterbi_LSIG = vitdec(LSIG_Deinterleave,trellis,1,'trunc','hard');
% Error_BCC_LSIG = sum(Viterbi_LSIG ~= L_SIG);
% 
% m = zeros(12,1);
% for i=1:1:12
%     m(i) = (2^(i-1))*Viterbi_LSIG(i+5);
% end
% DeL_LENGTH = (4*(sum(m)+3)/3)+20;
% 
% 
% %% VHT Signal B Decoding Processing
% 
% DeVHT_SIGB = Estimate_Data(577:640);
% for i = 1:1:64
%    if (DeVHT_SIGB(i) <=0)
%        DeVHT_SIGB(i) = 0;
%    else
%        DeVHT_SIGB(i) = 1;
%    end
% end
% 
% De_SIGB = zeros(52,1);
% for i = 1:1:7
%     De_SIGB(i) = DeVHT_SIGB(i+4);
%     De_SIGB(i+45) = DeVHT_SIGB(i+54);
% end
% for i = 1:1:13
%     De_SIGB(i+7) = DeVHT_SIGB(i+12);
%     De_SIGB(i+32) = DeVHT_SIGB(i+40);
% end
% for i = 1:1:6
%     De_SIGB(i+20) = DeVHT_SIGB(i+26);
%     De_SIGB(i+26) = DeVHT_SIGB(i+33);
% end
% Error_SIGB = sum(Interleaver_VHTSIGB ~= Interleaver_VHTSIGB);
% SIGB_DeinterleaveTable = Deinterleave_Table(:,1)+1;
% 
% SIGB_Deinterleave = zeros(52,1);
% for i = 1:1:52
%     n = SIGB_DeinterleaveTable(i);
%     SIGB_Deinterleave(n,1) = De_SIGB(i);
% end
% Error_SIGB_Deinterleave = sum(VHT_SIG_B_BCCEncode~=SIGB_Deinterleave);
% 
% Viterbi_SIGB = vitdec(SIGB_Deinterleave,trellis,1,'trunc','hard');
% Error_BCC_SIGB = sum(Viterbi_SIGB ~= VHT_SIG_B);
% 
% m = zeros(17,1);
% for i=1:1:17
%     m(i) = (2^(i-1))*Viterbi_SIGB(i);
% end
% DeSIG_B_LENGTH = sum(m);
% Error_SIGB_length = sum(VHT_SIG_B_LENGTH ~= DeSIG_B_LENGTH);
% De_APEPLengthMCS = DeSIG_B_LENGTH * 4;
% 
% %% VHT Signal A Decoding Processing
% 
% DeVHT_SIGA = Estimate_Data(321:448);
% 
% for i = 1:1:64
%     if (DeVHT_SIGA(i) <=0)
%         DeVHT_SIGA(i) = 0;
%     else
%         DeVHT_SIGA(i) = 1;
%     end
%     if ((DeVHT_SIGA(i+64)) <=0)
%         DeVHT_SIGA(i+64) = 0;
%     else
%         DeVHT_SIGA(i+64) = 1;
%     end
% end
% 
% DeVHTSIGA = zeros(48,2);
% for n =1:1:2
%     for i = 1:1:5
%         DeVHTSIGA(i,n) = DeVHT_SIGA(i+6+64*(n-1));
%         DeVHTSIGA(i+43,n) = DeVHT_SIGA(i+54+64*(n-1));
%     end
%     for i = 1:1:13
%         DeVHTSIGA(i+5,n) = DeVHT_SIGA(i+12+64*(n-1));
%         DeVHTSIGA(i+30,n) = DeVHT_SIGA(i+40+64*(n-1));
%     end
%     for i = 1:1:6
%         DeVHTSIGA(i+24,n) = DeVHT_SIGA(i+33+64*(n-1));
%         DeVHTSIGA(i+18,n) = DeVHT_SIGA(i+26+64*(n-1));
%     end
% end
% 
% SIGA_DeinterleaveTable = Deinterleave_Table(1:48,2)+1;
% 
% SIGA_Deinterleave = zeros(48,2);
% for k = 1:1:2
%     for i = 1:1:48
%         n = SIGA_DeinterleaveTable(i);
%         SIGA_Deinterleave(n,k) = DeVHTSIGA(i,k);
%     end
% end
% 
% SIGA_Deinterleave_map = [SIGA_Deinterleave(:,1); SIGA_Deinterleave(:,2)];
% Error_SIGA_Deinterleave = sum(SIGA_Deinterleave_map~=VHT_SIG_A_BCCEncode);
% 
% Viterbi_SIGA = vitdec(SIGA_Deinterleave_map,trellis,1,'trunc','hard');
% Error_BCC_SIGA = sum(Viterbi_SIGA ~= VHT_SIG_A );
% 
% De_BW = 2*Viterbi_SIGA(1) + Viterbi_SIGA(2);
% if (De_BW == 0)
%     Wifi_BW = 20;
% elseif (De_BW == 1)
%     Wifi_BW = 40;
% elseif (De_BW == 2)
%     Wifi_BW = 80;
% elseif (De_BW == 3)
%     Wifi_BW = 160;
% end

