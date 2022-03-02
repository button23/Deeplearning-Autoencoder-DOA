classdef QPSKBitsGenerator < matlab.System
    %#codegen
    % Generates the bits for each frame
    
    %   Copyright 2012-2017 The MathWorks, Inc.
    
    properties (Nontunable)
        ScramblerBase = 2;
        ScramblerPolynomial = [1 1 1 0 1];
        ScramblerInitialConditions = [0 0 0 0];
        NumberOfMessage = 10;
        MessageLength = 16;
        MessageBits = [];
    end
    
    properties (Access=private)
        pHeader
        pScrambler
        pSigSrc
    end
    
    properties (Access=private, Nontunable)
        pBarkerCode = [+1 +1 +1 +1 +1 -1 -1 +1 +1 -1 +1 -1 +1]; % Bipolar Barker Code
    end
    
    methods
        function obj = QPSKBitsGenerator(varargin)
            setProperties(obj,nargin,varargin{:});
        end
    end
    
    methods (Access=protected)
        function setupImpl(obj, ~)
            % Generate unipolar Barker Code and duplicate it as header
            ubc = ((obj.pBarkerCode + 1) / 2)';
            temp = (repmat(ubc,1,2))';
            obj.pHeader = temp(:);
            
            % Initialize scrambler system object
            obj.pScrambler = comm.Scrambler( ...
                obj.ScramblerBase, ...
                obj.ScramblerPolynomial, ...
                obj.ScramblerInitialConditions);
            
            % Initialize signal source
            obj.pSigSrc = dsp.SignalSource(obj.MessageBits, ...
                'SamplesPerFrame', obj.MessageLength * 7 * obj.NumberOfMessage, ...
                'SignalEndAction', 'Cyclic repetition');
            
        end
        
        function [y, msgBin] = stepImpl(obj)
            
            % Generate message binaries from signal source.
            msgBin = obj.pSigSrc();
            
            % Scramble the data
            scrambledMsg = obj.pScrambler(msgBin);
            
            % Append the scrambled bit sequence to the header
            y = [obj.pHeader ; scrambledMsg];
            
        end
        
        function resetImpl(obj)
            reset(obj.pScrambler);
            reset(obj.pSigSrc);
        end
        
        function releaseImpl(obj)
            release(obj.pScrambler);
            release(obj.pSigSrc);
        end
        
        function N = getNumInputsImpl(~)
            N = 0;
        end
        
        function N = getNumOutputsImpl(~)
            N = 2;
        end
    end
end

