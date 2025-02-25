% Mark 10:27 Jesus looked at them and said, "With man it is impossible,
% but not with God. For with God all things are possible."
%
% process_radar_sat1_data.m
% This program perform strip map SAR processing using a sample set of data
% from a RADARSAT-1 satellite. The Range Doppler Algorithm is used to process
% the data. The steps involved in this processing were taken from the book 
% "Digital Processing of Synthetic Aperture Radar Data" by Ian G. Cumming and 
% Frank H. Wong. ISBN-13: 978-1-58053-058-3. 
%
% The RDA will be performed using the following steps:
% 
% Range Compression -> Azimuth FFT -> RCMC -> Azimuth Compression
% -> Azimuth IFFT and Look Summation->Compressed Data
% 
% Each step is broken down into sub-steps which contain another set of
% steps specific to it. This program is written in a modular way to allow for
% easy reuse of certain parts of the code if needed. 
% 
% Range & Azimuth FFT: In the digital processor implementation, it does not
% matter whether the Range Fourier Transform or Azimuth Fourier Transform
% is performed first. In either case, the two-dimensional signal spectrum
% can be accurately represented by the 2D Frequency domain signal. Equation
% (5.24) on page 176 "Digital Processing of Synthetic Aperture Radar
% Data" book.

% Notes:
% microsecond = 1e-6, millisecond = 1e-3
% 1 MHz = 1e6,1 GHz = 1e9, 1KHz = 1e3
% 
% The term "range cell" refers to samples in the range or "fast
% time" direction. 
% 
% The term "range line" refers to samples in the azimuth or "slow time"
% direction.

%% Add paths to procesing scripts
clear,   home,   format compact
addpath(pwd + "\matlab-sar\sar-data-processor\sar-algorithms\range-doppler-algorithm\");


%% Perform pre-processing of data using script provided by book ISBN-13: 978-1-58053-058-3
specify_parameters

%% Init radarsat-1 data set params

% These params are from extract_dat.m and are specific to the data set a
% so they can not be changed.
radarsat1.length_replica  = length_replica; % Total length (I&Q) of replica record
radarsat1.tot_nrg_cells   = tot_Nrg_cells;  % Total number of range cells per line
radarsat1.tot_nrg_lines   = tot_Nrg_lines;  % Total number of range lines (records)
radarsat1.first_rng_line  = first_rg_line;  % first range line for this parameter set
radarsat1.first_rng_cell  = first_rg_cell;  % first range cell for this parameter set
radarsat1.first_replica   = first_replica;  % First record that contains the replica 
radarsat1.prf_hz          = PRF;            % Pulse Reputation Frequency (Hz)
radarsat1.fs_hz           = Fr;             % Radar sampling rate (Hz)
radarsat1.f0_hz           = f0;             % Radar center frequency (Hz)
radarsat1.c_mps           = c;              % Speed of light (m/s)
radarsat1.r0_m            = R0;             % Slant range of first radar sample (m)
radarsat1.nrepl           = Nrepl;          % No. of valid samples in the replica
radarsat1.kr_hz_s         = Kr;             % FM rate of radar pulse (Hz/s) (i.e. chirp rate)
radarsat1.tr_s            = Tr;             % Chirp duration (s)
radarsat1.vr_mps          = Vr;             % Effective Radar velocity (m/s)
radarsat1.n_rng_cells_fft = Nrg_cells;
radarsat1.n_rng_lines_fft = Nrg_lines;

% now clear cd params since we've saved them to a structure
clear_cd_params

%% Setup script 
% set to 1 to see plots, set to 0 for no plots
plot_real_sig = 0;
plot_raw_sig = 0;
plot_rng_comp_sig = 0;
plot_rng_dopp_sig = 0;
disp_baseband_dopp = 0;

% set to 1 to peform IFFT during range compression step, 0 otherwise
ifft_after_range_comp = 1;

%% calcuate important params
% calculate sweep bandwith for pulse, chirp rate * chirp duration
sweep_bandwidth_hz = radarsat1.kr_hz_s * radarsat1.tr_s;

% Calculate range resolution
range_resolution = physconst('LightSpeed')/(2*sweep_bandwidth_hz);

%% now print important radar parameters for this data set
fprintf(1,'------------------------------------ \n');
fprintf(1,'Radar Parameters \n');
fprintf(1,'Radar Center Frequency  \n\t%2.2f GHz \n',   radarsat1.f0_hz/1e9);
fprintf(1,'Radar Sampling Rate     \n\t%2.2f MHz \n',   radarsat1.fs_hz/1e6);
fprintf(1,'Radar PRF               \n\t%2.2f KHz \n',   radarsat1.prf_hz/1e3);
fprintf(1,'Radar Bandwidth         \n\t%2.2f MHz\n',    sweep_bandwidth_hz/1e6);
fprintf(1,'Radar FM Rate           \n\t%2.2f Hz/us\n',  radarsat1.kr_hz_s/1e12);
fprintf(1,'Radar Chirps            \n\t%2.2f  \n',      radarsat1.tot_nrg_lines);
fprintf(1,'Radar Chirp Duration    \n\t%2.2f us \n',    radarsat1.tr_s*1e6);
fprintf(1,'Radar Range Resolution  \n\t%2.2f m\n',      range_resolution);
fprintf(1,'Radar Velocity          \n\t%2.2f m/s\n',    radarsat1.vr_mps);
fprintf(1,'Radar Slant Range to 1st Sample \n\t%2.2f Km \n',radarsat1.r0_m/1e3);
fprintf(1,'Radar Range Cells       \n\t%2.2f \n',       radarsat1.n_rng_cells_fft);
fprintf(1,'Radar Range Lines       \n\t%2.2f \n',       radarsat1.n_rng_lines_fft);


%% 1. Range Compression
ref_chirp_signal = range_chirp(fs_hz = radarsat1.fs_hz, ...
                               pulse_duration_s = radarsat1.tr_s, ...
                               sweep_bandwidth_hz = sweep_bandwidth_hz, ...
                               number_of_samples = radarsat1.n_rng_cells_fft, ...
                               causal=0);
if(plot_real_sig)
    clc;close all;
    figure
    subplot(2,1,1)
    plot(real(ref_chirp_signal),'b')
    hold on
    plot(imag(ref_chirp_signal),'r')
    hold off
    grid on
    axis([0 radarsat1.n_rng_cells_fft -1.2 1.2])
    title(sprintf(['Real and Imaginay Parts of Time Series \n' ...
        'Fs= \t%2.2f Mhz, BW= \t%2.2f MHz, Pulse Duration= \t%2.2f usec'],...
        radarsat1.fs_hz/1e6, sweep_bandwidth_hz/1e6, radarsat1.tr_s*1e6));
    xlabel('Time Index(samples)')
    ylabel('Amplitude')
    
    subplot(2,1,2)
    plot(real(transpose(data)),'b')
    hold on
    plot(imag(transpose(data)),'r')
    hold off
    grid on
    axis([0 radarsat1.n_rng_cells_fft*2 -1.2 1.2])
    xlabel('Time Index')
    ylabel('Amplitude')
    title('SAR Received Signal')

end

% perform range compression
rng_compressed_sig = rng_compress_signal(ref_sig = ref_chirp_signal, ...
                                         recv_sig = data, ...
                                         perform_ifft = ifft_after_range_comp, ...
                                         number_of_samples = radarsat1.n_rng_cells_fft);

% Visualize results
if(plot_raw_sig)
    figure
    subplot(2,1,1)
    imshow(imread("RADARSAT_1.PNG"));
    ylabel('Cross-Range Samples')
    xlabel('Range Samples')
    title('Image region extracted from RADARSAT-1 data')

    subplot(2,1,2)
    imagesc(real(data));
    axis image              % adjust displayed aspect ratio
    legend
    colorbar;
    ylabel('Cross-Range Samples')
    xlabel('Range Samples')
    title('Real portion of I/Q data samples of extracted image region')

end

if(plot_rng_comp_sig)
    figure
    subplot(2,1,1)
    %Plot result to verify matched filter result
    plot(abs(fftshift(transpose(rng_compressed_sig))));
    axis([0 radarsat1.n_rng_cells_fft 0 80e3])
    colormap('gray')        % adjust the contrast
    grid on
    hold on
    xlabel('Samples')
    ylabel('Amplitude')
    title('Matched Filter Output')
    
    % display mangitude of range compressed data
    subplot(2,1,2)
    imagesc(abs(rng_compressed_sig));
    colormap('gray')        % adjust the contrast
    xlabel('Range bin');
    ylabel('Azimuth bin');
    title('Real Range Compressed Signal')
    
end

%% 2. Azimuth FFT
disp '2.Azimuth FFT'
% Azimuth FFT:
% An Azimuth FFT transforms the data into the range Doppler domain, where
% the Doppler Centroid estimation and most subsequent operations are
% performed. This is done by performing a FFT on every column of data.
% Recall for this data set rows = PRIs(Range lines), cols = Samples(Rage
% bins) so across PRI's

rng_dpplr_sig = fftshift(fft(rng_compressed_sig),1);
% Display iamgesc() of data
if(plot_rng_dopp_sig)
    figure(4)
    imagesc(abs(rng_dpplr_sig));
    axis image              % adjust displayed aspect ratio
    colormap('gray')        % adjust the contrast
    xlabel('Range bin');
    ylabel('Azimuth bin');
    title('Real Range Compressed Signal with Azimuth FFT')
end

%% 2.1 Doppler centroid estimation

% The doppler centroid is the "average doppler shift" or "center doppler frequency",
% formulated as f_nc = -Ka*n_c, where f_nc = doppler centroid, 
% Ka = Azimuth FM rate, and n_c = beam center crossing time relative to time
% of closest approach. This parameter must be known for further processing
% of SAR data.

% Define paraemeters
% The unaliased doppler centroid is given by: f_nc = f'nc + Mamb*Fa, where
f_nc = 0;       % f_nc  = Absolute doppler centroid frequency(Hz), 
f_prime_nc = 0; % f'n_c = Fracational PRF part of Baseband dopp freq (Hz),
Mamb = 0;       % Mamb  = Doppler Ambiguity #, integer multiple of Fa
Fa = radarsat1.prf_hz;         % Fa    = Pulse repetition frequency(Hz)

% In order to calculate the absolute doppler centroid we must first estimate
% 2 unknown parameters: f'n_c and Mamb. The fracational PRF part f'n_c sets 
% the azimuth matched filter center frequency and throwaway region. The ambiguity
% number indicates the position along the frequency axis of the absolute doppler
% centroid. The absolute centroid f_nc is used in RCM.
%
% Pg. 237: The RCM correction curve should be taken from azimuth
% frequencies lying within (+-)Fa/2 of the absolute doppler centroid, so that
% the correct curve is used, and the correction discontinuity is as far as 
% possible from the center of the main target energy. In order to estimate 
% the absolute azimuth frequency we must know the correct doppler centroid.
% Therefore we must calculate the absolute doppler centroid. 
%
% Doppler Centroid Estimation steps
%
% 1. Estimate the baseband component. There are 2 methods of doing
% this: Spectral fit method or Phase Increment method
%
% Perform f'n_c estimate using Phase Increment method.
%
% Doppler centroid estimation using Average phase increment method:
% The Doppler Centroid (DC) algorithms estimate the center frequency of the
% Doppler spectrum of the data, related to the azimuth beam center. The DC
% locates the azimuth signal energy in the azimuth frequency domain. It is 
% needed so that the signal energy in the Doppler spectrum can be correctly 
% captured by the azimuth compression filter, providing the best 
% signal-to-noise ratio and azimuth resolution.
% The fractional PRF f'nc sets the azimuth matched filter center frequency
% and throwaway region. 

%% Break up data in sub blocks and process 1 sub  block at a time

% TODO SDD: My need to create an outer loop that loads 1 block at at time
% nad processes data. 

Nrowsg = 3;    
Ncolsg = 3;    
Nspect = Nrowsg*Ncolsg;              % Total number of spectra to calculate
Nrglpb = floor( radarsat1.n_rng_cells_fft/Nspect );  % No. of range cells per spectra  
r1 = 1 + (1-1)*Nrglpb;   r2 = r1 + Nrglpb - 1;

data_for_est_baseband_dopp = zeros(size(data,1),length([r1:r2]), Nspect);

% Loop over each range line and extract a # of range cells
for range_line_idx = 1 : Nspect 
    % Calculate indexes of range cells that will be used to find the az
    % power spectra to be averaged( r1 = range cell 1, r2 = range cell 2
    % there are Nrglpb range cells per spectra
    r1 = 1 + (range_line_idx-1)*Nrglpb;   r2 = r1 + Nrglpb - 1;

    % Save range lines. Go from r1 to r2 which is a total of 227 range cells.
    data_for_est_baseband_dopp(:,:,range_line_idx) = data(:,r1:r2);
end 


% calculate the baseband dopp freq for each range line 1st then average them together. 
sub_block1_size = size(data_for_est_baseband_dopp(:,:,1),2);

% Loop over each sub-block and perform doppler centroid processing
f_prime_nc = zeros(sub_block1_size,1).';
accc_angle = zeros(sub_block1_size,1).';
fnc = zeros(1,Nspect);
for sub_blk_idx = 1:Nspect
    % extract sub block data
    sub_block_data = test_data_for_est_baseband_dopp(:,:,sub_blk_idx); % TODO SDD; modified this line to use "block 1" I got from compute_azim_spectra.m

    % Loop over each range sample
    for rng_sample_idx = 1:sub_block1_size

        % Extract range line from the sub-block
        range_line = sub_block_data(:, rng_sample_idx).';

        % Calculate baseband dopp for this range line
        [f_prime_nc(rng_sample_idx), accc_angle(rng_sample_idx)] = calc_centroid_baseband_madsens(slow_time_data=range_line,...
                                                    prf=Fa,...
                                                    plot="false");

        % check for wrapping of phase angle between azimuth looks
        % Calculate the phase difference
        phase_diff = diff(accc_angle);
        
        % Find indices where phase jump exceeds pi
        wrap_indices = find(abs(phase_diff) > pi, 1); % This returns the indices where
                                                   % the phase jumps exceed π, indicating wrapping points.

        if(~isempty(wrap_indices))
            warning("Phase Wrapping occured")
        else
            disp("No Phase Wrapping detected between ACCC angles")
        end

        if(disp_baseband_dopp)
            % find the mean of the sub-block to compare against spectral fit method
            % used in the compute_azim_spectra.m file. 
            fprintf('Spectral fit baseband dopp est for range cells: %d - %d: \n \t%2.2f \n', ...
                first_rg_cell, ...
                first_rg_cell + sub_block1_size-1, ...
                baseband_dopp_cent_solution(1,sub_block1_size));
            fprintf('Madsens baseband dopp est for range cells: %d - %d:\n  \t%2.2f \n', ...
                first_rg_cell, ...
                first_rg_cell + sub_block1_size-1, ...
                mean(f_prime_nc));
        end
    end
    % calculate baseband doppler centroid for this azimuth block
    baseband_dopp_centroid = mean(f_prime_nc);

    [dopp_amb, ~] = calc_dopp_ambiguity(Data=sub_block_data, ...
    F0_hz=radarsat1.f0_hz, ...
    Prf_hz=radarsat1.prf_hz,...
    F_prime_nc=f_prime_nc, ...
    Fm_rate_hz=radarsat1.kr_hz_s, ...
    Chirp_s=radarsat1.tr_s, ...
    Accc_angles=accc_angle, ...
    Plot_results=true);
    fprintf('Doppler Ambiguity Estimate w/ WDA: \t %2.2f \n', dopp_amb);
    
    % 3. Estimate the absolute doppler centroid for each range line in SAR data
    % calculate absolute dopp centroid 
    fnc(sub_blk_idx) = baseband_dopp_centroid; % + (dopp_amb * radarsat1.prf_hz);
    fprintf('Doppler Centroid Estimate Sub-Block %d: \t %2.2f \n',sub_blk_idx,fnc(sub_blk_idx));
end

% now loop over each range cell in the sub block and compute the baseband
% dopp frequency for that range line. This is effectively passing in a set 
% of azimuth samples per range cell from our data.
% test_baseband_dopp = test_data_for_est_baseband_dopp(:,:,1);
% % calculate the baseband dopp freq for each range line 1st then average them together. 
% sub_block1_size = size(test_data_for_est_baseband_dopp(:,:,1),2);
% f_prime_nc = zeros(sub_block1_size,1).';
% accc_angle = zeros(sub_block1_size,1).';
% for rng_sample_idx = 1:sub_block1_size
% 
%     % grab range line from the sub-block
%     range_line = sub_block_data(:, rng_sample_idx).';
% 
%     % calculate baseband dopp for this range line
%     [f_prime_nc(rng_sample_idx),accc_angle(rng_sample_idx)] = calc_centroid_baseband_madsens(slow_time_data=range_line,...
%                                                 prf=Fa,...
%                                                 plot="false");
% end

% check for wrapping of phase angle between azimuth looks
% Calculate the phase difference
% phase_diff = diff(accc_angle(:,1));
% 
% % Find indices where phase jump exceeds pi
% wrap_indices = find(abs(phase_diff) > pi); % This returns the indices where
%                                            % the phase jumps exceed π, indicating wrapping points.
% 
% if(~isempty(wrap_indices))
%     warning("Phase Wrapping occured")
% else
%     disp("No Phase Wrapping detected between ACCC angles")
% end
% 
% if(disp_baseband_dopp)
%     % find the mean of the sub-block to compare against spectral fit method
%     % used in the compute_azim_spectra.m file. 
%     fprintf('Spectral fit baseband dopp est for range cells: %d - %d: \n \t%2.2f \n', ...
%         first_rg_cell, ...
%         first_rg_cell + sub_block1_size-1, ...
%         baseband_dopp_cent_solution(1,1));
%     fprintf('Madsens baseband dopp est for range cells: %d - %d:\n  \t%2.2f \n', ...
%         first_rg_cell, ...
%         first_rg_cell + sub_block1_size-1, ...
%         mean(f_prime_nc(1,:)));
% end
% 
% % calculate baseband doppler centroid for this azimuth block
% baseband_dopp_centroid = mean(f_prime_nc);
% 
% % 2. Estimate the doppler ambiguity using a Doppler Ambiguity
% % Resolver(DAR). There are 2 approaches to doing this: Magnitude-based
% % methods or Phase-based methods our solution performs the doppler ambiguity
% % estimate(Mamb) using a Phase-based method.
% 
% num_az_s = size(sub_block_data,1);
% 
% % DAR needs the following
% % Average Cross Correlation Coefficient @ each range bin - ACCC
% % Pulse Repitition Frequency - PRF
% % Radar Center Frequency - F0
% 
% [dopp_amb, ~] = calc_dopp_ambiguity(Data=sub_block_data, ...
% F0_hz=radarsat1.f0_hz, ...
% Prf_hz=radarsat1.prf_hz,...
% F_prime_nc=f_prime_nc, ...
% Fm_rate_hz=radarsat1.kr_hz_s, ...
% Chirp_s=radarsat1.tr_s, ...
% Accc_angles=accc_angle, ...
% Plot_results=true);
% 
% fprintf('Doppler Ambiguity Estimate w/ WDA: \t %2.2f \n', dopp_amb);
% 
% % 3. Estimate the absolute doppler centroid for each range line in SAR data
% % calculate absolute dopp centroid 
% fnc = baseband_dopp_centroid + (dopp_amb * radarsat1.prf_hz);
% fprintf('Doppler Centroid Estimate Blocks 1-3: \t %2.2f \n',fnc(1:3));


%% important calculations
% Range of closest approach - R0
% r0_m = calc_rng_of_closest_approach(slant_range_m=data_freq(:,pulse_idx), ...
%                                     slow_time_secs=azimuth_time_sec, ...
%                                     vr_mps=vr_mps);
% % beam center crossing time - nc
% beam_cen_crss_time_secs = calc_beam_cent_crss_time(r0_m=r0_m, ...
%                                                    theta_deg=squint_angle_deg, ...
%                                                    vel_mps=vr_mps);
% 
% % slant range w/ respect to beam center crossing time - R(nc)
% slant_range_m = calc_slant_range(time=beam_cen_crss_time_secs, ...
%                                  type='beam center');