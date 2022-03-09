function est = vhtSingleStreamChannelEstimate(x,cfgVHT)
%vhtSingleStreamChannelEstimate Single stream channel estimates at pilot
%subcarrier location using VHT-LTF field
%
%   EST = vhtSingleStreamChannelEstimate(X,CFGVHT) returns the estimated
%   channel at pilot subcarrier locations for each VHT-LTF symbol, assuming
%   one space-time stream at the transmitter.
%   
%   EST is a complex Nsp-by-Nsym-by-Nr array characterizing the estimated
%   channel at pilot subcarrier locations for each symbol, where Nsp is the
%   number of occupied pilot subcarriers, Nsym is the number of VHT-LTF
%   symbols, and Nr is the number of receive antennas.
%
%   X is a complex Nst-by-Nsym-by-Nr array containing demodulated VHT-LTF
%   OFDM symbols. Nst is the number of occupied subcarriers.
%
%   CFGVHT is a packet format configuration object of type <a href="matlab:help('wlanVHTConfig')">wlanVHTConfig</a>. 

%   Copyright 2018 The MathWorks, Inc.

%#codegen

% Validate symbol type
validateattributes(x, {'single','double'}, {'3d'}, mfilename, 'VHT-LTF OFDM symbol(s)');

% validation object
validateattributes(cfgVHT, {'wlanVHTConfig'}, {'scalar'}, mfilename, 'VHT format configuration object');
    
chanBW = cfgVHT.ChannelBandwidth;
numSC = size(x,1);
numRxAnts = size(x,3);
numSTSTotal = sum(cfgVHT.NumSpaceTimeStreams);

% Return an empty if empty symbols
if isempty(x)
    est = zeros(numSC,numSTSTotal,numRxAnts);
    return;
end

% Get OFDM configuration
[cfgOFDM,~,pilotInd] = wlan.internal.wlanGetOFDMConfig(chanBW, 'Long', 'VHT', numSTSTotal);

% Channel estimate for single stream pilots
[ltf,Pheltf,numLTF] = wlan.internal.vhtltfSequence(chanBW, numSTSTotal);
R = Pheltf(1,1:numLTF); % R matrix changes pilot polarity per symbol

% Number of symbols in the input must be equal to numLTF
numSym = size(x,2);
coder.internal.errorIf(numSym~=numLTF, 'wlan:vhtSingleStreamChannelEstimate:InvalidOFDMSymbol', numSym, numLTF);

% Estimate the channel at pilot subcarriers accounting for polarity
if coder.target('MATLAB')
    est = x(pilotInd,:,:)./(ltf(cfgOFDM.PilotIndices).*R);
else % Codegen path
     est = coder.nullcopy(complex(zeros(numel(pilotInd),numLTF,numRxAnts)));
     for n=1:numRxAnts
         for m=1:numLTF
            est(:,m,n) = x(pilotInd,m,n)./(ltf(cfgOFDM.PilotIndices).*R(m));
         end
     end
end

end
