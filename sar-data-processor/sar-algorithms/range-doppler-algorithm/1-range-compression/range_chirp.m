function [chirp_signal,time_vector] = range_chirp(rangeChirp)
% range_chirp This function synthesize a linear FM chirp signal. 
% It is assumed to be a positive chirp. A chirp signal is a sinusoid whose 
% freuqency changes linearly with time. A formula for such a signal can be defined by
% creating a complex exponential with quadratic phase.
%
% Fs - Signal sample rate, specified as a positive scalar. Units are Hertz.
% 
% PulseDuration-  Length of each pulse (in seconds) as a positive scalar. The value must 
% satisfy PulseWidth <= 1./PRF.
%
% SweepBandwidth - FM sweep bandwidth in Hz
%
arguments
    rangeChirp.fs_hz              (1,1) {mustBePositive} = 0
    rangeChirp.pulse_duration_s   (1,1) {mustBePositive} = 0
    rangeChirp.sweep_bandwidth_hz (1,1) {mustBePositive} = 0
    rangeChirp.number_of_samples  (1,1) {mustBePositive} = 0
    rangeChirp.causal             (1,1) {mustBeInRange(rangeChirp.causal,0,1)} = 0
end
fs_hz             = rangeChirp.fs_hz;
pulse_duration    = rangeChirp.pulse_duration_s;
sweep_bandwidth   = rangeChirp.sweep_bandwidth_hz;
number_of_samples = rangeChirp.number_of_samples;
causal            = rangeChirp.causal;

if ~exist('causal','var')
     causal = 1;
end

if sweep_bandwidth/2 > fs_hz
    disp('Warning. Aliasing will be produced since BW/2 > fs');
end
    f0_hz = -sweep_bandwidth/2;
    f1_hz = sweep_bandwidth/2;
    if(causal)
        % Causal signal should be the same shape. Our starting point would be
        % zero instead of NumberOfSamples/2
        m =  0:number_of_samples - 1;

        m_squared = ((m) .* (m)) / 2;

        phi = ((f1_hz - f0_hz)/fs_hz) * (m_squared/number_of_samples) + (f0_hz/fs_hz) * (m);        

        chirp_signal = exp(1i*2*pi*phi);
        
        time_vector = 0:1/fs_hz:pulse_duration - (1/fs_hz);
    else
        % Divide by Fs to normalize frequency, spectrum goes to -Fs/2 to +Fs/2                                                            
        % Multiply by zero because we are non-causal 
        m =  -number_of_samples/2:number_of_samples/2 - 1;
        m_squared = ((m).*(m))/2;
        phi=((f1_hz - f0_hz)/fs_hz)* (m_squared/number_of_samples) + (f0_hz/fs_hz) * (m) *0;
        chirp_signal = exp(1i*2*pi*phi);
        time_vector = 0:1/fs_hz:pulse_duration - (1/fs_hz);
    end    

end

