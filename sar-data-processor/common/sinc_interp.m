% Mark 10:27 Jesus looked at them and said, "With man it is impossible,
% but not with God. For with God all things are possible."
function [reconstructed_signal] = sinc_interp(sincInterp)
% This function performs a sinc interpolation on a signal that meets the
% following conditions:
%   - The signal is bandlimited; that is, it's highest frequency is finite.
%     Measurements of any physical system are bandlimited
%
%   - The sampling satisfies Nyquist sampling rate. For a real signal, the
%     sampling rate must twice was high as the signals highest frequency
%     componenet. For a complex signal, the sampling rate must be higher than
%     the signals bandwidth.
%
% Inputs:
%   Ts_s             - sampling period in seconds
%   time_axis_s      - original time axis in seconds
%   interp_time_axis - interpolated time axis in seconds
%
% Outputs:
%   reconstructed_signal - The result of the matched filter operation (vector)
%
arguments
    sincInterp.sampled_signal     {mustBeNonempty} = 0
    sincInterp.Ts_s               {mustBeNonempty} = 0
    sincInterp.time_axis_s        {mustBeNonempty} = 0
    sincInterp.interp_time_axis_s {mustBeNonempty} = 0
end
t_interp_s = sincInterp.interp_time_axis_s;
x   = sincInterp.sampled_signal;
Ts  = sincInterp.Ts_s;
t_s = sincInterp.time_axis_s;


% For each discrete sample in the signal, generate a a new signal weighted 
% by a sinc function 
% Sinc Interpolation
reconstructed_signal = zeros(size(t_s)); % Initialize interpolated signal
for k = 1:length(t_interp_s)
    % Calculate interpolated value at each new time point and
    % sum the contributions from each original sample point across the entire
    % sampling period to create the reconstructed signal
    reconstructed_signal(k) = sum(x .* sinc((t_interp_s(k) - t_s) / Ts));
end

end