function [unambiguous_range] = calc_unamb_rng(unambigR)
%CALC_UNAMB_RNG Summary of this function goes here
%   Detailed explanation goes here
arguments
    unambigR.Prf_hz (1,1) {mustBeNonnegative}  = 0 
end
c = physconst('LightSpeed');
unambiguous_range = c./(unambigR.Prf_hz*2);
end

