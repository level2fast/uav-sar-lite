classdef RadarSat1Sar < Radar
    %RADARSAT1RADAR Class that describes RADARSAT1 parameters. This class
    % defines properties and methods specific to the RADARSAT1 radar
    % system.
    properties(Access = public)
        % Waveform properties
        Fnc_hz      (1,1) {mustBeNumeric}  = 0 % absolute doppler centroid TODO: CALC
        Frac_fnc_hz (1,1) {mustBeNumeric}  = 0 % fractional PRF part of doppler centroid TODO: CALC
        Kr_hz_s     (1,1) {mustBeNumeric}  = 0 % FM rate of transmitted pulse chirp, pos for up chirp, neg for down chirp
        Km_hz_s     (1,1) {mustBePositive} = 0 % Range FM rate of point target signal in the range doppler domain(combining Kr and Ksrc) TODO: CALC
        Ka_hz_s     (1,1) {mustBePositive} = 0 % Azimuth FM rate of point target signal TODO: CALC
        Ksrc_hz_s   (1,1) {mustBeNumeric}  = 0 % FM rate of secondary compression filter
        La_m        (1,1) {mustBePositive} = 0 % Antenna aperture in azimuth direction
        Ls_m        (1,1) {mustBePositive} = 0 % Synthetic aperture 
        R0_m        (1,1) {mustBeNumeric}  = 0 % Slant range of closest approach TODO: CALC
        Vr_mps      (1,1) {mustBeNumeric}  = 0 % Effective radar velocity, meters per second
        Vg_mps      (1,1) {mustBeNumeric}  = 0 % Antenna beam footprint, meters per second  TODO: CALC
        Vr_ref_mps  (1,1) {mustBeNumeric}  = 0 % Effective radar velocity at reference range, meters per second
        Vrel_mps    (1,1) {mustBeNumeric}  = 0 % Relative velocity between radar and target, meters per second  TODO: CALC
        Vs_mps      (1,1) {mustBeNumeric}  = 0 % Radar platform velocity, meters per second  TODO: CALC
        Ta_s        (1,1) {mustBeNumeric}  = 0 % Azimuth exposure time  TODO: CALC
        Theta_cc_rad(1,1) {mustBeNumeric}  = 0 % Range/azimuth cross coupling phase term TODO: CALC
        Theta_i_rad (1,1) {mustBeNumeric}  = 0 % incidence angle TODO: CALC
        Theta_n_rad (1,1) {mustBeNumeric}  = 0 % Off-nadir angle TODO: CALC
        Theta_r_rad (1,1) {mustBeNumeric}  = 0 % General squint angle in rectilinear geometry, measured from zero doppler TODO: CALC
        Theta_rc_rad(1,1) {mustBeNumeric}  = 0 % Squint angle of radar beam center in rectilinear geometry, measured from zero doppler TODO: CALC
        Waveform    (1,:) {mustBeNonempty} = 0 % Radar waveform used for transmit and processing
    end

    % methods
    %     function outputArg = calcGroundResolution(object)
    %         %METHOD1 Summary of this method goes here
    %         %   Detailed explanation goes here
    %         outputArg = object;
    %     end
    %     function outputArg = calcGroundBasis(object)
    %         %METHOD1 Summary of this method goes here
    %         %   Detailed explanation goes here
    %         outputArg = object;
    %     end
    %     function outputArg = calcIntegrationTime(object)
    %         %METHOD1 Summary of this method goes here
    %         %   Detailed explanation goes here
    %         outputArg = object;
    %     end
    %     function outputArg = calcSlantRangeToGroundRange(object)
    %         %METHOD1 Summary of this method goes here
    %         %   Detailed explanation goes here
    %         outputArg = object;
    %     end
    %     function outputArg = calcCrossRangeResolutionMeters(object)
    %         %METHOD1 Summary of this method goes here
    %         %   Detailed explanation goes here
    %         outputArg = object;
    %     end
    %     function outputArg = calcCrossRangeResolutionSeconds(object)
    %         %METHOD1 Summary of this method goes here
    %         %   Detailed explanation goes here
    %         outputArg = object;
    %     end
    %     function range_resolution = calcRangeResolution(object)
    %         %METHOD1 Summary of this method goes here
    %         %   Detailed explanation goes here
    %         range_resolution = object;
    %     end
    % end

end




