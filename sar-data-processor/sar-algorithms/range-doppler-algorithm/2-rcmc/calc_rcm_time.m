function [rcm_time] = calc_rcm_time(lambda_m,R0,fn_hz,Vr_mps)
%CALC_RCM_SAMPLES Summary of this function goes here
%   Detailed explanation goes here
rcm_meters = lambda_m^2 * R0 * fn_hz / 8 *Vr_mps;
rcm_time = (2 * rcm_meters) /physconst('LightSpeed');
end

