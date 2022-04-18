

clear;close all;clc;
numOFDM = 1;
Nfft = 1024 ;
% Nfft = 10^5;
Ts = 100/1024*exp(-6);
M = 256;
k = log2(M);
cpLength = 30;
% cpLength = 30:-4:10;
delay = 10;
EbNo = 0:30;
% EbNo = 20;


PAPR = zeros(21,Nfft);
ber = zeros(length(cpLength),length(EbNo));
numError = zeros(size(EbNo)); 

h = [1 0.8];
A = sqrt(1/170);

llr = zeros(k,1);
demodulatedData = zeros(k,Nfft*numOFDM/k);

modTx = zeros(Nfft,numOFDM);




%%

for pp = 1 : length(cpLength)
    for bb = 1 : length(EbNo)
        % Convert Eb/No to SNR
        snrdb = EbNo(bb) + 10*log10(k);
        % Noise variance calculation for unity average signal power
        noiseVar = 10.^(-snrdb/10);
        %% QPSK modulation
        dataIn = randi([0 1],Nfft*k,numOFDM); 
        compareIn = reshape(dataIn,8,Nfft*numOFDM);
        modTx = qammod(dataIn,M,'InputType','bit','UnitAveragePower',true);   
%         modTx2 = qammod(dataIn,M,'InputType','bit');    
%         modTxTest = modTx2 / sqrt(170); 
        


%         %% IFFT
%         txIfft = ifft(modTx);
%         papr = abs(txIfft);
%         sumPapr = sum(papr,2)/numOFDM;
%         normPapr = sumPapr/max(sumPapr);
% 
%         %% Add CP
%         txCP = [(txIfft(end-cpLength(pp)+1:end,:));txIfft];

        %% Parallel to stream
%         stream = reshape(txCP,(Nfft+cpLength(pp))*numOFDM,1);
        stream = reshape(modTx,(Nfft)*numOFDM,1);

        

        %% Multipath 
%         path1 = stream;
%         path2 = [zeros(delay,1);stream(1:end-10)];
%         mulPath = 0.8 * path1 + 0.4 * path2;
        
        test_data = conv(stream,h);
        test_data = test_data(1:Nfft);
        
        h_zero = [h,zeros(1,Nfft-length(h))];
        hfft = fft(h_zero);
        Hk = (abs(hfft));

        %% AWGN 
%         rxSig = awgn(mulPath,snrdb,'measured');
        rxSig = awgn(test_data,snrdb,'measured');


        %% Stream to Parallel
%         paraData = reshape(rxSig,Nfft+cpLength(pp),numOFDM);
        paraData = reshape(rxSig,Nfft,numOFDM);


%         %% Remove CP
%         rxRemCpData = paraData(cpLength(pp)+1:end,:);

%         %% FFT
%         rxFFT = fft(rxRemCpData);

        %% QPSK demodulation (Soft value) 
%         a =1;
        for oo = 1 : numOFDM
        for ll = 1: Nfft 
            rr = real(paraData(ll,oo));
            im = imag(paraData(ll,oo));
            llr(1) = rr; 
            llr(2) = - abs(rr) + 8 * A;
            llr(3) = - abs( abs(rr) - 8 * A ) + 4 * A;
            llr(4) = - abs( abs( abs(rr) - 8 * A ) - 4 * A) + 2 * A;
            
            llr(5) = -im; 
            llr(6) = - abs(im) + 8 * A;
            llr(7) = - abs( abs(im) - 8 * A ) + 4 * A;
            llr(8) = - abs( abs( abs(im) - 8 * A ) - 4 * A) + 2 * A;
            
            llr = llr * Hk(ll)^2 ; 

            
            demodulatedData(:,ll+(oo-1)*Nfft) = llr > 0;
        end
           
        end
         numError(bb) = biterr(demodulatedData,compareIn);
         ccc = demodulatedData - compareIn;
%         %% QPSK demodulation
%         rxData = qamdemod(rxFFT,M,'OutputType','bit', ...
%                     'UnitAveragePower',true,'NoiseVariance',noiseVar);
%                 
%         %% Calculate bit errors            
%          numError(bb) = biterr(rxData,dataIn);
         
    end
     ber(pp,:) = numError/(Nfft * numOFDM);
end
    
    %% Plot the result
    semilogy(EbNo,ber(1,:),'-oc')
%     hold on
%     semilogy(EbNo,ber(2,:),'-*b')
%     semilogy(EbNo,ber(3,:),'-+m')
%     semilogy(EbNo,ber(4,:),'-sk')
%     semilogy(EbNo,ber(5,:),'-dg')
%     semilogy(EbNo,ber(6,:),'-hr')

%     legend('CP = 10','CP = 8','CP = 6','CP = 4','CP = 2','CP = 0')
%     legend('CP = 30','CP = 26','CP = 22','CP = 18','CP = 14','CP = 10')
    grid
    xlabel('Eb/No(dB)')
    ylabel('Bit Error Rate')
% 
%     figure
%     plot(normPapr)
%     xlabel('delay')
%     ylabel('ber')
%     xlim([500 1000])
%     
%     xlabel('Sample')
%     ylabel('Amplitude of IFFT Result(Normalized)')