function [range_vec,doppler_vec,slow_time_vec] = create_rng_dopp_vectors(radar)
%CREATEGRIDS Summary of this function goes here
%   Detailed explanation goes here
arguments
    radar.Fs       (1,1) {mustBeNonnegative}  = 0  % sampling rate
    radar.Prf      (1,1) {mustBeNonnegative}  = 0 % pulse repitition frequency
    radar.NumPulse (1,1) {mustBeNonnegative}  = 0 % number of pulses in a CPI
end
n_samp_pri = round(radar.Fs/radar.Prf); % number of samples for 1 pulse repitition interval (PRI)
delta_r = physconst('LightSpeed')/(2*radar.Fs); % represents spacing between range grid samples
range_vec = (0:n_samp_pri-1) *delta_r; % number of samples within a PRI spaced by delta_r

% create a vector with values that are spaced by a fraction or % of our PRF
% Since pulses are viewed from the perspective of time we could say that the 
% slow time vector indices represent a fraction of the time of each pulse 
% we transmit.Recall that the time step across slow time is PRI=1/PRF. 
slow_time_vec = (0:(radar.NumPulse-1))/radar.Prf; 


delta_f = radar.Prf/radar.NumPulse;
doppler_vec = -radar.Prf/2:delta_f:((radar.Prf/2)-delta_f); % doppler grid samples, always between (-prf/2):(prf/2)
end

