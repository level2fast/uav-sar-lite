% Mark 10:27 Jesus looked at them and said, "With man it is impossible,
% but not with God. For with God all things are possible."
function [reconstructed_signal] = sinc_interp2(samp_sig, orig_fs, new_time_axis)
%sinc_interp Summary of this function goes here
% This function performs a sinc interpolation on a signal that meets the
% following conditions
%   -The signal is bandlimited; that is, it's highest frequency is finite.
%   Mesurements of any physical system are bandlimited
%   -The sampling satisfies Nyquist sampling rate. For a real signal, the
%   sampling rate must twice was high as the signals highest frequency
%   componenet. For a complex signal, the sampling rate must be higher than
%   the signals bandwidth.
%
% sampSig:
%   Sampled signal
%
% origFs:
%   Sampling rate of the orignal signal
%
% newTimeAxis:
%   Axis we are interpolating into

% New signal
t = new_time_axis'; % time axis

% Original sampled signal
Ts = 1/orig_fs;  % sampling period 
xn = samp_sig;   % sampled signal

sinc_train = zeros(length(t), length(xn));
sig_len = length(xn);
nind = -length(xn)/2 + length(xn)/2 +1;
win = kaiser(length(sig_len),2.5);
for n = -floor(length(xn)/2):floor(length(xn)/2)
    if(nind < sig_len + 1)
        sinc_train(:, nind) = xn(nind) * sinc((t - n * Ts) / Ts) .* win;
    end
    nind = nind + 1;
end

xr = sum(sinc_train, 2); % sum(sincTrain, 2) is the interpolated/reconstructed 
reconstructed_signal = xr;
end
