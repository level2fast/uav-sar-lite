
function [dopp_amb, remainder] = calc_dopp_ambiguity(doppAmbRes)
%CALC_DOPP_AMBIGUITY Estimates the doppler ambiguity of an input signal
%   The function of the doppler ambiguity resolver (DAR) is to estimate the
%   correct ambiguity number, Mamb. There are 2 types of DAR's broadly broken
%   up into 2 categories. Magnitude based algorithms and phased based algorithms.
%   
%%  Phased Based DAR
%   Wavelength Diversity Algorithm (WDA)
%       Operates on range frequency data. Uses the slant range as a function
%       of the beam center crossing time nc, to calculate the derivative of
%       the ACCC with respect to range frequency(ft), which is then used to
%       get slope of the angle of the ACCC which is then used to find the 
%       absolute doppler frequency.
%
%   Inputs
%      Data            - Range look vector. The second dimension must contain 
%                        2 looks when using MLCC or MLBF.
%      F_prime_nc_hz   - Baseband doppler frequency in hertz.
%      Prf_hz          - Pulse repitition frequency.
%      F0_hz           - Radar carrier frequency or center frequency.
%      Type            - Doppler ambiguity estimation method: Wavelength Diversity
%                        Algorithm, Multilook Cross Correlatoin Algorithm,
%                        Multilook Beat Frequency Algorithm.
%  
%  Outputs
%      dopp_amb      - doppler ambiguity 
arguments
    doppAmbRes.Data       (:,:) = 0
    doppAmbRes.F_prime_nc {mustBeNonempty} = 0
    doppAmbRes.Prf_hz     {mustBeNonempty} = 0
    doppAmbRes.F0_hz      {mustBeNonempty} = 0
    doppAmbRes.Fm_rate_hz {mustBeNonempty} = 0
    doppAmbRes.Chirp_s    {mustBeNonempty} = 0
    doppAmbRes.Accc_angles {mustBeNonempty} = 0
    doppAmbRes.Type string {mustBeMember(doppAmbRes.Type,{'WDA','MLCC','MLBF'})} = 'WDA'
    doppAmbRes.Plot_results logical = 0
end
type           = doppAmbRes.Type;
data           = doppAmbRes.Data;
f_prime_nc_vec = doppAmbRes.F_prime_nc;
prf_hz         = doppAmbRes.Prf_hz;
f0_hz          = doppAmbRes.F0_hz;
accc_angles    = doppAmbRes.Accc_angles.';
Kr             = doppAmbRes.Fm_rate_hz;
Tr             = doppAmbRes.Chirp_s;
plot_results   = doppAmbRes.Plot_results;
if(isequal(type,'MLCC') | isequal(type ,'MLBF'))
    looks = size(data);
    if(looks(2) ~= 2)
        error("Input Data must contain 2 looks");
    end
end
% get size of range dim
len = size(data,2);
% calculate the slope of the ACCC angle vs. range frequency using linear
% regression(least squares approx).
ones_col = (ones(len,1));
range_freq_sweep = Kr*Tr/2;
rng_start_idx = -(range_freq_sweep);
rng_end_idx = range_freq_sweep;
range_axis = linspace(rng_start_idx, rng_end_idx, len).';  % in Hz
data_freq = [ones_col range_axis]; % use ones to get y-int when perfomring linear regression
accc_slope =  data_freq \ accc_angles;
% Extract slope and intercept
y_intercept = accc_slope(1);
slope = real(accc_slope(2));
if(y_intercept > 0)
    warning("ACCC vs. Range Freq y-int %d", y_intercept)
end
% calculate absolute dopp centroid 
f_prime_nc = mean(f_prime_nc_vec);
% estimate the absolute dopp freq using slope of ACCC
fnc_hz = (prf_hz/(2 * pi)) * f0_hz * slope;
% calculate the doppler ambiguity to within an integer value of the PRF
dopp_amb = round((fnc_hz - f_prime_nc) / prf_hz);
% calculate the remainder
% this value is a meausure of the accuracy of the absolute dopp freq, as an
% example of an accuracy criterion the remainder should be less than 33% of the PRF.
remainder = dopp_amb * prf_hz - (fnc_hz - f_prime_nc);
if(remainder < 0.33)
    warning("Remainder should be less than 33% of a PRF")
end
if(plot_results)
    % Plot the angle vector against the frequency vector
    disp(accc_slope(1));
    disp(accc_slope(2));
    remainder_perc = (remainder/prf_hz) * 100;
    % Define slope and intercept
    m = slope;       % Slope
    b = y_intercept;       % Y-intercept
    
    % Calculate y values based on slope and intercept
    y = m * range_axis + b;
    
    figure;
    plot(range_axis/1e6, accc_angles, 'b-', 'LineWidth', 1.5);
    hold on
    grid on
    plot(range_axis/1e6,y)
    xlim([min(range_axis/1e6) max(range_axis/1e6)])
    xlabel('Range Frequency (MHz)');
    ylabel('ACCC Angle (radians )');
    t = sprintf(['ACCC Angle vs. Frequency \n ' ...
        'Slope of Line Fit = %f mrad/Mhz, \t ' ...
        'Intercept =   %2.2f, \t ' ...
        'Dopp Amb Est =  %2.2f, \t' ...
        ' Remainder =  %2.2f %% of PRF'], ...
        (m), ...
        (b), ...
        (dopp_amb), ...
        (remainder_perc));
    title(t);
end
end