function [doppler_amb] = calc_centroid_dopp_amb(avgPhaseChng)
%CALC_CENTROID_DOPP_AMB Summary of this function goes here
%   Detailed explanation goes here
arguments
    avgPhaseChng.data (:,:) {mustBeNonempty}
    avgPhaseChng.prf  (1,1) {mustBeNonempty}
    avgPhaseChng.plot (1,1) {mustBeInRange(avgPhaseChng.plot,0,1)} = 0
end
doppler_amb = avgPhaseChng.data;

end

