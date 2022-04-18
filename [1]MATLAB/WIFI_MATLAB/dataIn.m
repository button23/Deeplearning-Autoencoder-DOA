% 
% n =1000;
% indata_real  = randi([-n n],64,1);
% indata_imaginary = randi([-n n],64,1);
% indata_fft = indata_real +1i*indata_imaginary;

% indata_fft_wifi = fft_data_real*1000 +1i*fft_data_imaginary*1000;
% real_scale_part = int16(real(indata_fft_wifi));
% real_scale_imaginary = int16(imag(indata_fft_wifi));
% indata_fft_wifi = real_scale_part + 1i*real_scale_imaginary;
% 


% indata_fft_wifi = []
% fft_data= ((fft( indata_fft)));
% data_fft_shift = dsp_test_fftshift;
% % data_fft_shift = fftshift(fft_data,1);
% fft_data_scale0 =int16(data_fft_shift * 8);
% fft_data_scale1 =int16(data_fft_shift * 4);
% fft_data_scale2 =int16(data_fft_shift * 2);
% fft_data_scale3 =int16(data_fft_shift );
% fft_data_scale4 =int16(data_fft_shift/2 );
% fft_data_scale5 =int16(data_fft_shift/4 );
% fft_data_scale6 =int16(data_fft_shift/8 );


%% scrambler
% intscram =[ 1 0 1 1 1 0 1];
regis = ones(7,1);
% datain = rand([0 1],127,1);
sequence = zeros(127,1);
for ss = 1:127
    tem = xor(regis(1),regis(4)); 
    regis(1:end-1) = regis(2:end);
    regis(end) = tem;
    sequence(ss) = tem;
end