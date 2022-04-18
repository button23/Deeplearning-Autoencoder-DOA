%%
%
%  Read data from collected RF data by USRP
%  Creator: ZZF
%

%%  Basic Parameters
function output = IPS_dataLoad(Path, indexPoint, sampleFileName)
% To complete the path (find the folder name)
numFolder = num2str(indexPoint);
% Number of collections for each point
intePath = strcat(Path, numFolder, sampleFileName);
f=fopen(intePath);
dataa = textscan(f, '%d16');
fclose(f);
output = dataa{1};
end