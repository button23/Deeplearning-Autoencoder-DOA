function [nest,sigest] = vhtNoiseEstimate(x,chanEstSSPilots,cfgVHT,varargin)
%vhtNoiseEstimate Estimate noise power using VHT data field pilots
%
%   NEST = vhtNoiseEstimate(X,CHANESTSSPILOTS,CFGVHT) estimates the mean
%   noise power in watts using the demodulated pilot symbols in the VHT
%   data field and single-stream channel estimates at pilot subcarriers.
%   The noise estimate is averaged over the number of symbols and receive
%   antennas.
%
%   X is the received time-domain VHT Data field signal, specified as an
%   Ns-by-Nr matrix of real or complex values. Ns represents the number of
%   time-domain samples in the VHT Data field and Nr represents the number
%   of receive antennas.
%
%   CHANESTSSPILOTS is a complex Nsp-by-Nltf-by-Nr array containing the
%   channel gains at pilot subcarrier locations for each symbol, assuming
%   one space-time stream at the transmitter. Nsp is the number of pilots
%   subcarriers and Nltf is the number of VHT-LTF symbols.
%
%   CFGVHT is the format configuration object of type <a
%   href="matlab:help('wlanVHTConfig')">wlanVHTConfig</a>.
%   
%   NEST = vhtNoiseEstimate(...,SYMOFFSET) specifies the sampling offset as
%   a fraction of the cyclic prefix (CP) length for every OFDM symbol, as a
%   double precision, real scalar between 0 and 1, inclusive. The OFDM
%   demodulation is performed based on Nfft samples following the offset
%   position, where Nfft denotes the FFT length. The default value of this
%   property is 0.75, which means the offset is three quarters of the CP
%   length.
%
%   [NEST,SIGEST] = vhtNoiseEstimate(...) additionally returns an estimate
%   of the signal power.

%   Copyright 2018 The MathWorks, Inc.

%#codegen

% Validate inputs
validateattributes(x, {'double'}, {'2d','finite'}, mfilename, 'VHT-Data field input'); 
validateattributes(chanEstSSPilots, {'double'}, {'3d','finite'}, mfilename, 'Single stream pilots channel estimate input'); 
validateattributes(cfgVHT, {'wlanVHTConfig'}, {'scalar'}, mfilename, 'VHT format configuration object');

if nargin == 4
    symOffset = varargin{1};
    % Validate symbol offset
    validateattributes(symOffset, {'double'}, {'real','scalar','>=',0,'<=',1}, mfilename, 'OFDM sampling offset');
else
    symOffset = 0.75;
end

numSTS = sum(cfgVHT.NumSpaceTimeStreams);
cfgOFDM = wlan.internal.wlanGetOFDMConfig(cfgVHT.ChannelBandwidth, cfgVHT.GuardInterval, 'VHT', numSTS);

% Ensure at least 1 OFDM symbol in the input
numSamples = size(x,1);
minMultipleLength = cfgOFDM.FFTLength+cfgOFDM.CyclicPrefixLength;
numOFDMSym = floor(numSamples/minMultipleLength);
minInputLength = minMultipleLength;
coder.internal.errorIf(numSamples<minInputLength, 'wlan:shared:ShortDataInput', minInputLength);

% OFDM demodulate
minInputLength = minMultipleLength*numOFDMSym;
[~,ofdmDemodPilots] = wlan.internal.wlanOFDMDemodulate(x(1:minInputLength,:), cfgOFDM, symOffset);

% Get reference pilots, from Eqn 22-95, IEEE Std 802.11ac-2013
% Offset by 4 to allow for L-SIG, VHT-SIG-A, VHT-SIG-B pilot symbols
n = (0:numOFDMSym-1).';
z = 4; 
% Set the number of space time streams to 1 since the pilots are same
% across all spatial streams
refPilots = wlan.internal.vhtPilots(n, z, cfgVHT.ChannelBandwidth, 1);

% Estimate CPE and phase correct symbols
% Average single-stream pilot estimates over symbols (2nd dimension)
chanEstSSPilotsAvg = mean(chanEstSSPilots,2);
[cpe,estRxPilots] = wlan.internal.commonPhaseErrorEstimate(ofdmDemodPilots, chanEstSSPilotsAvg, refPilots);
ofdmPilotsData = wlan.internal.commonPhaseErrorCorrect(ofdmDemodPilots, cpe);

% Estimate noise
pilotError = estRxPilots-ofdmPilotsData;
nest = mean(real(pilotError(:).*conj(pilotError(:))));

if nargout>1
    % Get power of channel estimate at pilot locations
   sigest =  mean(estRxPilots(:).*conj(estRxPilots(:)));
end

end
