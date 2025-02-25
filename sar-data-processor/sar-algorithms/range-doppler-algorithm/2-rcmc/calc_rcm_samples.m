function [rcm_samples] = calc_rcm_samples(lambda_m,R0,fn_hz,Vr_mps,fs_hz)
%CALC_RCM_SAMPLES Summary of this function goes here
%   Detailed explanation goes here
rcm_meters = lambda_m^2 * R0 * fn_hz / 8 *Vr_mps;
rcm_samples = (rcm_meters * fs_hz * 2) /physconst('LightSpeed');
end

