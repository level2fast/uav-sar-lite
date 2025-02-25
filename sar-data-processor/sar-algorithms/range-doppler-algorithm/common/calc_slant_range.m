function [outputArg1] = calc_slant_range(slantRange)
%CALC_SLANT_RANGE Calculates slant range in meters from radar to target
% 
%   Inputs
%      time_secs      -
%      type - 
%  
%  Outputs
%      beam_center_crossing_time_secs - 
arguments
    slantRange.time_secs = 0
    slantRange.type = 'beam center'
end
type = slantRange.type;
time = slantRange.time;

if(isequal(type,'beam center'))
    % calculate slant range as function of the beam center crossing
    % time
end

end

