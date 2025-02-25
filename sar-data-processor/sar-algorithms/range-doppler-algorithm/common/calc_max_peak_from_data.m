function rca_rng = calc_max_peak_from_data(range_data, num_rng_bins)
    % calculate_RCA_from_data - Calculates the range of closest approach (RCA) for each target
    %
    % Syntax:
    % RCA_ranges = calculate_RCA_from_data(range_data, range_resolution)
    %
    % Inputs:
    %   range_data - 2D matrix of range-compressed data (range bins x azimuth positions)
    %   range_resolution - Range resolution of the radar (meters per bin)
    %
    % Outputs:
    %   RCA_ranges - Array containing RCA for each target in meters
    
    % Find the number of range bins and azimuth positions
    [num_range_bins, num_azimuth_positions] = size(range_data(1:num_rng_bins, :));
    
    % Initialize an empty array to store RCA for each target
    rca_ranges = [];
    rca_peaks  = [];
    
    % Loop through each azimuth position to find the range bin with max intensity
    for az = 1:num_azimuth_positions
        % Find the peak in the current azimuth (across range bins)
        [max_intensity, range_bin] = max(abs(range_data(1:num_range_bins, az)));
        
        % Calculate the actual range of closest approach
        rca = range_bin;
        
        % Store RCA in meters
        rca_ranges = [rca_ranges; rca];
        rca_peaks = [rca_peaks; max_intensity];
    end

    [peak_intensity, peak_idx] = max(rca_peaks);
 
    % the range of closest approach should be the target with the
    % smallest distance among all possible ranges between the radar
    % platform and each target and the highest peak magnitude
    % rca_rng = min(rca_ranges); % Returns single RCA value per target across azimuths
    rca_rng = rca_ranges(peak_idx);   
    
    % Display RCA results
    fprintf('Range sample index of With Greatest Intensity: %.2f meters\n', rca_rng);
end
