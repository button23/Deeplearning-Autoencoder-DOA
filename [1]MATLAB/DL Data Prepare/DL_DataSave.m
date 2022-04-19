% Save the generated train data and test data to the specified path
% DATE: 2022.01.28 by ZZF
%
% Input:
%   dataType:       Training data or test data
%   operSys         The operating system used by the computer
%   snr             The signal to noise ratio (for naming data folders)
%   vec_all:        The vector data vectorized from noisy SCM
%   vec_all_ori:    The vector data vectorized from noiseless SCM
%   angle_doa:      The random angle (# sample x # source)
%
% Output:           None
%
function DL_DataSave(dataType,operSys,snr,vec_all,vec_all_ori,angle_doa)
snr = num2str(snr);
snr = append(snr,'dB');
if isequal(operSys, 'WINDOWS')
    %% WINDOWS
    rootPath = 'C:\Users\HYPC300\OneDrive - 한양대학교\GitHub\Deeplearning-Autoencoder-DOA\data';
elseif isequal(operSys, 'UBUNTU') %% UBUNTU
    rootPath = '/home/hymc/[0]Github/data';
else %% MAC
    rootPath = '/Users/button/Deeplearning-Autoencoder-DOA/data'; 
end

dataName = append(dataType,'_data');
oriDataName = append(dataType,'_origin_data');
labelName = append(dataType,'_label');

% create train or test folder 
basePath = fullfile(rootPath,date,snr,dataType); % fullfile() can be replaced by append(). data is a keyword in MATLAB.
if ~isfolder(basePath)                % add snr for differentiating the folder names
    mkdir(basePath);
end

%NOTE: The last string in the path is the .mat file name,
% but the variable inside the file might be different.
dataPath   = fullfile(basePath,dataName);
originPath = fullfile(basePath,oriDataName);
labelPath  = fullfile(basePath,labelName);

%% Saving data to the specific path
if isequal(dataType, 'train')
    % Save training data
    train_data = vec_all.'; % tranpose is to make the data conform to numsample by (number antenna element^2(noisy))
    train_origin_data = vec_all_ori.';
    train_label = angle_doa;
else
    % Save test data
    test_data = vec_all.';
    test_origin_data = vec_all_ori.';
    test_label = angle_doa;
end

% NOTE: The second argument is the name of the variable (must be rounded by ' ')
save(dataPath, dataName) % Save noisy data.
save(originPath, oriDataName) % Save noiseless data.
save(labelPath, labelName) % Save label (DOA angle)
end

