reference=textread('C:\Users\user\Desktop\wlan toolbox\m_files\referenceData.txt','%n','delimiter',' ');
rawdata=textread('C:\Users\user\Desktop\wlan toolbox\m_files\802.11ac.txt','0x%s','delimiter',',');

deci_data = hex2dec(rawdata(1:1024));
data_bin = dec2bin(deci_data);
data_bin_re = reshape(data_bin',8,1024*4);
data_bin_re =data_bin_re';
data_bin_no_sign = data_bin_re(:,2:end);
data_bin_int8 = bin2dec(data_bin_no_sign);



 deci_reference = hex2dec(reference);
reference_bin = dec2bin(deci_reference);
reference_bin_lr = fliplr(reference_bin);
reference_bin_reshape = reshape(reference_bin_lr',1,64*32);
reference_bin_reshape=reference_bin_reshape';
reference_bin_reshape_int8 = int8(reference_bin_reshape)-48;


% data = randi([0 1],372,1);
% data = ones(372,1);
data = [1;0;1;0;1;0;1;0;1;0;1;0;1;0;1;0];
data1= [data;data;data;data;data;data;1;0;1;0];
% data2 = [data1;data1;data1;data1;data1;data1;data1];
% data3 =[data2;data;data;data;data;1;0;1;0;1;0;1;0;];

trellis = poly2trellis(7,[133 171]);
rate = 1/2;
codedData = convenc(data3(1:372),trellis);
data_scale = -127*codedData+63;
decodedData = vitdec(codedData,trellis,34,'trunc','soft');



data_bcc_encode = wlanBCCEncode(reference_bin_reshape_int8,1/2);

 test = randi([0 1],60,1);
 test_real = 1-2*test;
 BM = zeros(30,2);
 a =1;
 for ii = 1:2:60
    BM(a,1)=test_real(ii)+test_real(ii+1);
    BM(a,2)=test_real(ii)-test_real(ii+1);
    a=a+1;
 end
  test_bcc_de = wlanBCCDecode(test,1/2);
% 
%  
% data_bcc_decode = wlanBCCDecode(data_bin_int8,1/2);


