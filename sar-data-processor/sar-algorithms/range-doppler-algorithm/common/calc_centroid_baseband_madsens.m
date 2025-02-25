function [baseband_dopp_centroid, accc_angle] = calc_centroid_baseband_madsens(doppCentroid)
%CALC_CENTROID_MADSENS_ALG: This function calculates the baseband doppler centroid.
%
%   The baseband doppler centroid sets the azimuth matched filter frequency 
%   and throwaway region. In this method, the phase of the complex radar
%   data is utilized to estimate the fractional part of the Doppler
%   Centroid. This is possible because the expected value of the azimuth
%   phase increment equals the phase increment at the center of the targets
%   exposure.
%
%   slow_time_data - range compressed data
%   prf  - pulse repetition frequency
%   plot - plot baseband doppler centroid
arguments
    doppCentroid.slow_time_data (1,:) {mustBeVector} = []
    doppCentroid.prf  (1,1) {mustBeNonempty}
    doppCentroid.plot  {mustBeMember(doppCentroid.plot,{'true','false'})} = 'false'
end
data = doppCentroid.slow_time_data;
prf_hz  = doppCentroid.prf;
plot_centroid = doppCentroid.plot;

% calcluate average cross correlation coefficient
accc = calc_average_cross_correlation_coeff("slow_time_data", data);

% Calculate the angle of the Average Cross Correlation Coeffeicient
accc_angle = angle(accc);

% f'nc or baseband doppler centroid is the fractional part of the absolute 
% doppler centroid equation which is fnc = f'nc + PRF * Mamb. 
f_prime_nc_hz = (accc_angle * prf_hz) / (2 * pi); % baseband doppler centroid [Hz]

if(plot_centroid == "true")

    all_sig_inc_vec = zeros([1,data_len]);
    for idx = 1:(data_len - 1)
        all_sig_inc_vec(idx) = (conj(data(idx)) * data(idx+1));
    end

    % create vector that includes sum of all signal incrment vectors
    all_sig_inc_vec = [all_sig_inc_vec accc];

    % normalize vectors to have good scaling
    all_sig_inc_vec = normalize(all_sig_inc_vec,"norm",Inf);

    % get vector size
    all_sig_inc_vec_size = length(all_sig_inc_vec);
    
    figure
    % Loop over all vectors plotting the magnitude and phase of each one.
    % Should see a large vector that stands out from the rest. This vector
    % indicates the portion the PRF that can be used to calculate the
    % baseband doppler centroid
    for idx2 = 1:(all_sig_inc_vec_size)
 
        % Extract the angles and normalized magnitudes
        inc_vec_angle = angle(all_sig_inc_vec(idx2));
        inc_vec_magnitude = abs(all_sig_inc_vec(idx2));

        [x, y] = pol2cart(inc_vec_angle,inc_vec_magnitude);
        plot(x,y,"o")
        ax = gca;
        ax.XAxisLocation = 'origin';
        ax.YAxisLocation = 'origin';
        axis equal
        xlabel("Real signal icrement vector")
        ylabel("Imaginary signal icrement vector")
        grid on
        hold on
        line([0, x], [0, y]);
    end
    msg = sprintf("Phase Increments from a single range line:\n f'nc - \t%2.2f" + ...
        " or \n \t%2.2f * PRF",f_prime_nc_hz, accc_angle/(2*pi));
    title(msg)
    hold off
end

baseband_dopp_centroid = f_prime_nc_hz;
end

% For more info see:
%    Madsens Average Phase Increment method. T See page 511 of "Digital Processing
%    of Synthetic Aperture Radar Data" by Ian G. Cumming and Frank H. Wong. ISBN-13: 978-1-58053-058-3.