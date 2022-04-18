%% clear
% clear
% clc
close all
%% Generate sequence
% seq = randi([1,10],64,1);
% n = 0:1023;


%% Zero Padding
k= 32;
seq_zeroPad = zeros(64*k,1);
seq_zeroPad((2048-64)/2+1:(2048-64)/2+64) = seq;
% seq_zeroPad(1:64) = seq;
test_seq = repmat(seq,32,1);

%% ifft
ifft_seq_ = ifft(seq_zeroPad);
% ifft_seq_ = ifft(ifftshift(seq_zeroPad,1));
test_seq(64*10+1:end) = 0;
ifft_seq_rep = ifft(test_seq);


ifft_seq_reshape = reshape(ifft_seq_,k,numel(ifft_seq_)/k);
ifft_seq_64 = ifft(seq);

expp1 = exp(2*pi*1i*1*n/1024);
expp2 = exp(2*pi*1i*20*n/1024);
expp3 = exp(2*pi*1i*100*n/1024);



%% Plot
%  k = 1;
% figure(1)
% stem(abs(ifft_seq_64))
% figure(2)
% stem(abs(ifft_seq_))
figure(3)
stem(abs(ifft_seq_rep))
% figure(4)
% stem(abs(ifft_seq_reshape(1,:)))
ratio = abs(ifft_seq_reshape(1,:).')./abs(ifft_seq_64);
%  figure(3)
%  subplot(3,1,1)
%  plot((expp1(1:k)))
%  axis([-1 1 -1 1])
%
%  subplot(3,1,2)
%  plot((expp2(1:k)))
%   axis([-1 1 -1 1])
%
%  subplot(3,1,3)
%  plot((expp3(1:k)))
%  axis([-1 1 -1 1])
%  legend("1","2","3")
