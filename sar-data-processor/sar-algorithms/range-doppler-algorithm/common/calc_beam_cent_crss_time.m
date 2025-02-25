function [beam_center_crossing_time_secs] = calc_beam_cent_crss_time(beamCenter)
%CALC_BEAM_CENT_CRSS_TIME Calculates beam center crossing time in seconds
%   This calculation is valid for scenarios that have low squint and moderate
%   aperture length. Like ERS or RADARSAT.
%       
%   Inputs
%      r0_m      - range of closest approach in meters
%      theta_deg - squint angle in degrees
%      vel_mps   - velocity in meters per second
%  
%  Outputs
%      beam_center_crossing_time_secs - time that the beam center crosses
%                                       the target in seconds
arguments
    beamCenter.r0_m      {mustBeNumeric} = 0;
    beamCenter.theta_deg {mustBeNumeric} = 3.9;
    beamCenter.vel_mps   {mustBeNumeric} = 7062;
end
r0_m      = beamCenter.r0_m;
theta_deg = beamCenter.theta_deg;
vel_mps   = beamCenter.vel_mps;

beam_center_crossing_time_secs = -(r0_m * tand(theta_deg))/vel_mps;
end

