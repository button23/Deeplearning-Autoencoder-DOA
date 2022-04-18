chNoise_ = comm.AWGNChannel('NoiseMethod','Signal to noise ratio (SNR)',...
    'SNR', -9);

% insignal = randi([0 1],100,1);
insignal_r = [1,2;3,4;5,6];
insign_i = [2,3;8,9;5,8];
vvv = insignal_r + j*insign_i;
%Send the input signal through the channel.
sssss = vvv(:,1);
outsignal = chNoise_(sssss);


rxSig_VHT = chNoise_(tgacChan(wholeWaveform_VHT));  %% wholeWaveform_scaling  wholeWavefor

rx_cc = (tgacChan(wholeWaveform_VHT));  %% wholeWaveform_scaling  wholeWavefor
rx_cc_f = fft(rx_cc);
trans = fft(wholeWaveform_VHT);
HHH = rx_cc_f./trans;

siv_withoutNoise = HHH .* trans;
siv_withoutNoise_ifft = ifft(siv_withoutNoise);

% ccccc = ones(1760,1);
% hhhh = tgacChan(ccccc);
% HH = fft(hhhh);
%  HH = hhhh;

sig_conv =conv(wholeWaveform_VHT,hhhh);
resss = sig_conv(1:1760);

noisee = rxSig_VHT - siv_withoutNoise_ifft;

noisee_s = sum((abs(noisee).^2));
snrrr2 = sum((abs(wholeWaveform_VHT)).^2)./noisee_s;

snrrr1 = 10*log10((snrrr2));

