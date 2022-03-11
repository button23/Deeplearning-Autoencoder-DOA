%% Save Path
savePath = 'C:\Users\HYPC300\Desktop\new\WIFI_bpsk_ips.txt';
saveBinPath = 'C:\Users\HYPC300\Desktop\new\WIFI_bpsk_ips.bin';

%% Concatenate real part and imaginary part
frame_sync_int16  = int16(wholeWaveform_scaling); 
real_data = fi(real(frame_sync_int16));
imag_data = fi(imag(frame_sync_int16));
cplx_data = bitconcat(real_data,imag_data);

%% Write data to bin file
% Text file open (the file needs to be created beforehand under following path)
f=fopen(savePath, 'r');
% Scan data from the opened file into the variable 'data' 
data=fscanf(f, '%u', [length(cplx_data),1]);
% Create the bin file under such path
f2=fopen(saveBinPath, 'w');
% Write data into bin file
fwrite(f2,data, 'uint32');
% Open the bin file and read the data from it to check if anything wrong
f2=fopen(saveBinPath, 'rb');
data2=fread(f2, [length(cplx_data), 1], 'uint32');

% Check if the written data is the same as the original data
isequal(cplx_data, data2)
