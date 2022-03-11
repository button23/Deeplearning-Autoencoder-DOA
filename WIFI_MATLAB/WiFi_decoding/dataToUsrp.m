
fftc_input_dsp = zeros(30720,1);
% #0 symbol: Scaling 2^20
fftc_input_dsp(1:2208) = subframe_data(1:2208) * 2^20;
% #1 symbol: Scaling 2^21
fftc_input_dsp(2209:4400) = subframe_data(2209:4400)  * 2^21;
% PDSCH symbol(Traffic data): Scaling 2^19
fftc_input_dsp(4401:end) = subframe_data(4401:end)  * 2^19;
% Convert to integer
fftc_input_dsp_int16  = int16(fftc_input_dsp);


eqGrid_scale_num = reshape(eqGrid_scale,numel(eqGrid_scale),1);
res_eqGrid_scale_16 = int16(eqGrid_scale_num);

%% wifi data from matlab to usrp
frame_sync_int16  = int16(wifi30(:, 2)); %% WIFI_Waveform  wholeWaveform_scaling  wifi30

real_data = fi(real(frame_sync_int16));
imag_data = fi(imag(frame_sync_int16));

cplx_data = bitconcat(real_data,imag_data);

% ??? ?? text file open
% f=fopen('C:\Users\HYPC300\Desktop\new\wifi_2_11.txt','r');
f=fopen('.\wifi_2_11_ant2.txt','r');

% Open? text file? data ??? ??
data=fscanf(f,'%u',[length(frame_sync_int16),1]);
% ???? ??? bin file open
% f2=fopen('C:\Users\HYPC300\Desktop\new\wifi_2_11.bin','w');
 f2=fopen('.\wifi_2_11_ant2.bin','w');
% data? ?? binary? write
fwrite(f2,data,'uint32');
% ??? ????? bin??? ?? ???? ??
% f2=fopen('C:\Users\HYPC300\Desktop\new\wifi_2_11.bin','rb');
 f2=fopen('.\wifi_2_11_ant2.bin','rb');
data2=fread(f2, [length(frame_sync_int16),1], 'uint32');
isequal(data, data2)



