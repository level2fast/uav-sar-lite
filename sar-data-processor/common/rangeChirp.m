% Mark 10:27 Jesus looked at them and said, "With man it is impossible,
% but not with God. For with God all things are possible."
function [chirpSignal,timeVector] = rangeChirp(Fs,PulseDuration,SweepBandwidth,NumberOfSamples,causal)
%LINEARFMCHIRP Summary of this function goes here
% This function synthesize a linearFMChirp signal. It is assumed to be
% a positive chirp.A chirp signal is a sinusoid whose freuqency changes
% linearly with time. A formula for such a signal can be defined by
% creating a complex exponential with quadratic phase.
% TODO SDD: Remember when using fft take the normalized fft of the waveform i.e.
% : waveform = waveform(:)/norm(waveform(:));
% Fs:
% Signal sample rate, specified as a positive scalar. Units are Hertz.
% 
% PulseDuration:
% Length of each pulse (in seconds) as a positive scalar. The value must 
% satisfy PulseWidth <= 1./PRF.
%
% SweepBandwidth:
% FM sweep bandwidth in Hz
%
% if ~exist('plot','var')
%     plot = 0;
% end
% if ~exist('useSamples','var')
%     useSamples = 0;
% end
if ~exist('causal','var')
     causal = 1;
end
 
if nargin<4
    error('Error. Not enough input arguments.')
end

if SweepBandwidth/2 > Fs
    disp('Warning. Aliasing will be produced since BW/2 > fs');
end
    f0 = -SweepBandwidth/2;
    f1 = SweepBandwidth/2;
    if(causal)
        % Causal signal should be the same shape. Our starting point would be
        % zero instead of NumberOfSamples/2
        M =  0:NumberOfSamples-1;
        M_squared = ((M).*(M))/2;
        phi=((f1 - f0)/Fs)* (M_squared/NumberOfSamples) + (f0/Fs) * (M);        
        chirpSignal = exp(1i*2*pi*phi);
        timeVector = 0:1/Fs:PulseDuration;
    else
        % Divide by Fs to normalize frequency, spectrum only goes to -Fs/2 to +Fs/2                                                            
        % Multiply by zero because we are non-causal 
        M =  -NumberOfSamples/2:NumberOfSamples/2;
        M_squared = ((M).*(M))/2;
        phi=((f1 - f0)/Fs)* (M_squared/NumberOfSamples) + (f0/Fs) * (M) *0;
        chirpSignal = exp(1i*2*pi*phi);
        timeVector = 0:1/Fs:PulseDuration;
    end    

end

