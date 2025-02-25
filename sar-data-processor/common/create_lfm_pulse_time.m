function [lfm_signal] = create_lfm_pulse_time(lfmPulse)
%CREATE_LFM_PULSE_SAMPLES Summary of this function goes here
%   This function synthesizes a linear FM pulse. An LFM pusle is a sinusoid
%   whose frequency changes linearly from some low value to a high value or
%   vice versa. The formula for such a signal can be represented as a
%   complex exponential with quadratic phase:
%   lfm_pulse = e^(j*s(t))
%       where
%           s(t) = 2*pi*alpha*t^2 + 2*pi*fc*t+phi
%       and
%       fs = signal sample rate. Specified as a positive scalar (Hz)
%       F0 = Start frequency (Hz)
%       F0 = End frequency (Hz)
%       T = Total duration of pulse
%       chirpUpDown = indicates a positive frequency slope (+1) or a negative frequency slope (-1)
arguments
    lfmPulse.Fs (1,1) {mustBeNonnegative}  = 0 
    lfmPulse.F0 (1,1) {mustBeNonnegative}  = 0 
    lfmPulse.F1 (1,1) {mustBeNonnegative}  = 0 
    lfmPulse.T  (1,1) {mustBeNonnegative}  = 0 
    lfmPulse.ChirpUpDown (1,1) {mustBeNonnegative}  = 0 
end
fs_hz = lfmPulse.Fs;
sweep_bandwidth = fs_hz * lfmPulse.T;
f1_hz = lfmPulse.F1;
f0_hz = lfmPulse.F0;
pulse_duration_s = lfmPulse.T;
phi = 0;
if(sweep_bandwidth/2 > fs_hz)
    disp('Warning: aliasing will be produced since BW/2 > fs')
end

if(lfmPulse.ChirpUpDown==-1)
    temp = lfmPulse.F0;
    lfmPulse.F0 = lfmPulse.F1;
    lfmPulse.f1=temp;
end

dt = 1/fs_hz;
time_vec_s = 0:dt:lfmPulse.T-(1/fs_hz);
k = (f1_hz - f0_hz) / pulse_duration_s;
freqz_hz = 2 * pi * ((k/2) .* time_vec_s + f0_hz) .* time_vec_s + phi;
theta = freqz_hz + phi;
lfm_signal = exp(1j*theta);
end

