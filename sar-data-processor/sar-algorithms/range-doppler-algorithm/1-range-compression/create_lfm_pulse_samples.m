% Mark 10:27 Jesus looked at them and said, "With man it is impossible,
% but not with God. For with God all things are possible."
function [lfm_signal] = create_lfm_pulse_samples(lfm)
%CREATE_LFM_PULSE_SAMPLES Summary of this function goes here
%   This function synthesizes a linear FM pulse. An LFM pusle is a sinusoid
%   whose frequency changes linearly from some low value to a high value or
%   vice versa. The formula for such a signal can be represented as a
%   complex exponential with quadratic phase:
%   lfm_pulse = e^(j*s(t))
%       where
%           s(t) = 2*pi*alpha*t^2 + 2*pi*fc*t+phi
%       and
%       fs              = signal sample rate. Specified as a positive scalar (Hz)
%       sweepBandwidth  = length of the frequncy ramp in (Hz)
%       numberOfSamples = Total amount of samples used
%       causal          = specifies domain of signal. causal(0), non-causal(1)
%       chirpUpDown     = indicates a positive frequency slope (+1) or a negative frequency slope (-1)
arguments
    lfm.Fs_hz            (1,1) {mustBeNonnegative}  = 0 
    lfm.SweepBandwidth_hz(1,1) {mustBeNonnegative}  = 0 
    lfm.NumberOfSamples  (1,1) {mustBeNonnegative}  = 0 
    lfm.ChirpUpDown      (1,1) {mustBeNonnegative}  = 0 
    lfm.Causal           (1,1) {mustBeNonnegative}  = 0 
end
phi = 0;

if(lfm.SweepBandwidth/2 > lfm.Fs)
    disp('Warning: aliasing will occur since BW/2 > fs')
end

if(lfm.ChirpUpDown == -1)
    f0_hz = (-lfm.SweepBandwidth/2) * -1;
    f1_hz =  (lfm.SweepBandwidth/2) * -1;
else
    f0_hz = (0);
    f1_hz = (lfm.SweepBandwidth);
end

if(lfm.Causal == 0)
    % total duration of pulse(secs)
    pulse_duration_s = lfm.NumberOfSamples/lfm.Fs;
    % our starting point would be zero instead of number samples/2
    time_vec_s = 0:1/lfm.Fs:pulse_duration_s-(1/lfm.Fs);
else % non-causal
    % N is our sample vector which counts the number of samples in this pulse.
    % Our starting point is numberOfSamples/2. The -1 accounts for MATLABS
    % indexing
    N = -lfm.NumberOfSamples/2:lfm.NumberOfSamples/2-1;

    % time vector is number of samples / sampling rate = secs/pulse
    time_vec_s = N/lfm.Fs;
    pulse_duration_s = lfm.NumberOfSamples/lfm.Fs;
end

k = (f1_hz - f0_hz) / pulse_duration_s;
freqz_hz = 2 * pi * ((k/2) .* time_vec_s + f0_hz) .* time_vec_s + phi;
theta = freqz_hz + phi;
lfm_signal = exp(1j*theta);
end

