% Estimate and compensate the phase offfset 
% DATE: 2021.12.09 by ZZF
%
% Input:
%   data:                     Received data from USRP
%   f_normalize:              Normalization flag  (1: normalize, 0: no normalize)

% Output:
%   pfo:                      Estimated phaseoffset

function [pfo]=phase_corr(data,f_normalize)
n_stream = length(data(1,:));
NormalizedData = data;
% Normalize the data to make the data of all the channel have the
% same amplitude.
if f_normalize
    amp = zeros(n_stream,1);
    for i = 1:n_stream
        amp(i) = max(abs(data(:,i)));
    end
    
    maxAmp = max(amp);
    for i = 1:n_stream
        NormalizedData(:,i) = maxAmp/amp(i)*data(:,i);
    end
    disp('Data being Normalized !')
    
end

% Estimate the phase difference between channels.
    freq = zeros(size(data));
    pfo = cell(n_stream,1);
    for i = 1:n_stream
        freq(:,i) = fft(NormalizedData(:,i));
    end
    
    for i = 1:n_stream-1
        % phase offset estimation
        ang = rad2deg(angle(max(freq(:,1))/max(freq(:,i+1))));
        pfo{i} = comm.PhaseFrequencyOffset('PhaseOffset',ang);
    end
    disp('Phaseoffset being estimated !')
end


