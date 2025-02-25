function [range_of_closest_approach] = calc_rng_of_closest_approach(rngCloseAppr)
%CALC_RNG_OF_CLOSEST_APPROACH Calculate the range of closest approach in meters.
% 
%   This function should find the smallest distance among all possible
%   ranges between the platform and each target when given a test 
%   data set with targets at varying ranges. The range of closest
%   approach acts as a reference point to help align phase
%   measurements for SAR in processing steps further down the pipeline. 
%
%   Range of closest approach refers to the shortest distance between the 
%   radar platform (e.g., an aircraft or satellite) and the target area on 
%   the ground. It corresponds to the point at which the platform is closest
%   to the target, typically occurring when the radar beam is pointed directly
%   perpendicular to the platform's flight path.
% 
%   The range of closest approach defines the point where the Doppler frequency
%   shift of the target is zero, which is essential for resolving the target's
%   location and for focusing the SAR image correctly. The range of closest 
%   approach in SAR (Synthetic Aperture Radar) processing is typically 
%   treated as a scalar quantity. It refers to the shortest distance between
%   the radar platform (aircraft, satellite, or UAV) and the target on the
%   ground during the SAR acquisition. This distance is essential in SAR 
%   geometry because it marks the point of minimum slant range, which helps
%   determine parameters like resolution and the range reference used for
%   focusing the image.
%
%   Inputs
%      slant_range_m  - distance from radar to target in meters 
%      vr_mps         - velocity of radar in meters per second
%      slow_time_secs - time vector that represents the azimuth postition along the flight path 
%  
%  Outputs
%      range_of_closest_approach - the line-of-sight distance between the 
%                                  radar and the ground target.
arguments
    rngCloseAppr.slant_range_m  (1,1) {mustBeScalarOrEmpty} = 0;
    rngCloseAppr.vr_mps         (1,1) {mustBeScalarOrEmpty} = 0;
    rngCloseAppr.slow_time_secs (:,:) {mustBeVector} = 0;
end
slant_range_m  = rngCloseAppr.slant_range_m;
vr_mps         = rngCloseAppr.vr_mps;
slow_time_secs = rngCloseAppr.slow_time_secs.';
vel_per_rng_bin = (vr_mps .* slow_time_secs) ./ slant_range_m;
range_of_closest_approach_vec = slant_range_m .* sqrt(1 - (vel_per_rng_bin).^2);
range_of_closest_approach = min(range_of_closest_approach_vec);
%range_of_closest_approach = sqrt(slant_range_m^2 - (vr_mps^2 .* slow_time_secs.^2));
end

