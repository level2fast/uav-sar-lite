 classdef RangeCellMigrationCorrectionTest < matlab.unittest.TestCase
    methods (Test)
        function testDataGen(testCase)
            % The first concept to understand when correcting for range
            % walk is how the range variation effects the time axis. Range
            % walk can be viewed across PRI's and one way to produce it in
            % a linear fashion in a simulation is by performing a range
            % frequency shift. A range slow time vector is needed to do
            % this. This test verifies that we've generated a range slow
            % time vector that displays a linear increase in time as a
            % result of a constant range and velocity factor.

            % Paramers
            fs_hz = 1e6;
            T = 10e-6;
            t = 0:1/fs_hz:T-1/fs_hz;
            v = 1000;
            R0 = 50;

            % simulate range walk
            range_slow_time_rw = R0 + v *t;
            range_slow_time    = R0 + t;

            % The following figure plots the range slow time vector with
            % and without a velocity change across 1 pulse. The plots
            % demonstrate the effect that a constant velocity has on the
            % slow time vector over the duratoin of a pulse.
            clf
            figure(100)
            plot(t*1e6, range_slow_time)
            hold on
            plot(t*1e6, range_slow_time_rw)
            title("Range Slow Time vs. Pulse Duration: R0 = 50m");
            xlabel("Pulse Duration (\mus)")
            ylabel("Range Slow Time (\mus)")
            lgd = legend("Range Slow Time Vec w/ Vel=0mps", "Range Slow Time Vec w/ Vel=1000mps");
            lgd.FontSize = 12;
            lgd.FontWeight = "bold";

            % now we check that range slow time is linearly increasing due
            % to motion of the target. Should see a millisecond increase in
            % range slow time vector due to effects of velocity
            range_slow_time_diff = range_slow_time_rw - range_slow_time;
            time_diff = diff(round(range_slow_time_diff,4)); % round to nearest millisecond
            actual = true;
            tolerance = 9e-4;
            expected = all((time_diff) >= tolerance);
            testCase.verifyEqual(actual, expected);
        end

        function testFrequencyShifting(testCase)
            import matlab.unittest.TestCase
            import matlab.unittest.constraints.IsEqualTo
            import matlab.unittest.constraints.AbsoluteTolerance
            import matlab.unittest.constraints.RelativeTolerance
            % The second concept to understand is the Fourier Transform
            % Time-Shift property which states: A shift of ‘𝒕𝟎’ in time domain
            % is equivalent to introducing a phase shift of –𝝎𝒕𝟎. But amplitude
            % remains same.This property is important because range walk
            % effects are essentialy a shift in the time of where the energy from
            % received signal shows up after being being received by the system.
            % The shift of the data in range across pulses can be simulated
            % using a phase shift applied to the range frequency portion of the
            % signal, but first we need to verify that we can perform frequency
            % shifting of a signal using the FT time-shift property. This test
            % performs that verification.

            % Simple script that demonstratse the fourier transform frequency shift
            % property which says that a shift s(t-t0) == S(omega) * S(omega * t0), where
            % omega = 2*pi*f

            % Parameters
            Fs_hz = 1e3 ;              % Sampling frequency (Hz)
            T = 1 / Fs_hz;             % Sampling interval (s)
            N = 1024;               % Number of points for FFT (power of 2 for efficiency)
            t = (0:N-1) * T;        % Time vector (0 to N-1 samples)
            % Frequency axis
            freq_axis = (0:N-1)*(Fs_hz/N);
            % freq_axis = freq_axis - Fs_hz*(freq_axis > Fs_hz/2); % Center frequencies at 0 Hz

            % Create an impulse response (e.g., a sinc function)
            fc = 50; % Cut-off frequency (Hz)
            impulse_response = sinc(2 * fc * (t - N*T/2)); % Centered sinc function

            % Compute the Fourier Transform of the impulse response
            freq_response = fft(impulse_response, N);

            % FWD Shift in the frequency domain
            freq_shift_fwd_hz = -100;       % Frequency shift (Hz)
            shift_exp_fwd = exp(1j * 2 * pi * freq_shift_fwd_hz * t); % Complex exponential
            shifted_freq_response_fwd = freq_response .* shift_exp_fwd; % Apply shift

            % Inverse Fourier Transform to get shifted impulse response
            shifted_impulse_response_fwd = ifft(shifted_freq_response_fwd, 'symmetric');

            % BACK Shift in the frequency domain
            freq_shift_back_hz = 100;       % Frequency shift (Hz)
            shift_exp_back = exp(1j * 2 * pi * freq_shift_back_hz * t); % Complex exponential
            shifted_freq_response_back = freq_response .* shift_exp_back; % Apply shift

            % Inverse Fourier Transform to get shifted impulse response
            shifted_impulse_response_back = ifft(shifted_freq_response_back, 'symmetric');

            % Check shift by comparing peak locations
            [~, peakIndexOriginal] = max(abs(impulse_response));
            [~, peakIndexShiftedFwd] = max(abs(shifted_impulse_response_fwd));
            [~, peakIndexShiftedBack] = max(abs(shifted_impulse_response_back));

            frequencyShiftFwd = freq_axis(peakIndexShiftedFwd) - freq_axis(peakIndexOriginal);
            disp(['Frequency shift Fwd: ', num2str(frequencyShiftFwd), ' Hz']);

            frequencyShiftBack = freq_axis(peakIndexShiftedBack) - freq_axis(peakIndexOriginal);
            disp(['Frequency shift Back: ', num2str(frequencyShiftBack), ' Hz']);

            actual = frequencyShiftFwd;
            expected = 100;
            tolerance = 1;
            % Verify with absolute tolerance
            testCase.verifyThat(actual, IsEqualTo(expected, 'Within', ...
                matlab.unittest.constraints.AbsoluteTolerance(tolerance)));

            % Plot original and shifted impulse responses in the time domain
            figure(3);
            subplot(3, 1, 1);
            plot(t, impulse_response, 'b-', 'LineWidth', 1.5);
            xlabel('Time (s)');
            ylabel('Amplitude');
            title('Original Impulse Response');
            grid on;

            subplot(3, 1, 2);
            plot(t, shifted_impulse_response_fwd, 'r-', 'LineWidth', 1.5);
            xlabel('Time (s)');
            ylabel('Amplitude');
            title(['Shifted Impulse Response (Frequency Shift = ' num2str(freq_shift_fwd_hz) ' Hz)']);
            grid on;

            subplot(3, 1, 3);
            plot(t, shifted_impulse_response_back, 'r-', 'LineWidth', 1.5);
            xlabel('Time (s)');
            ylabel('Amplitude');
            title(['Shifted Impulse Response (Frequency Shift = ' num2str(freq_shift_back_hz) ' Hz)']);
            grid on;
        end

        function testRangeMigrated2dDataGenTest(testCase)
            % Create a range migrated signal using a simple waveform. Range
            % migration should be clearly visible across entire CPI.

            % Paramters
            fs_hz = 40e6; % Sampling Frequency Hz
            fc_hz = 16e9; % Carrier Frequency Hz
            prf_hz = 60e3; % Carrier Frequency Hz
            pulse_width_s = 1/prf_hz; % Pulse width
            c_mps = physconst('LightSpeed'); % Speed of Light mps
            R0_m = 200;

            wvl = c_mps/fc_hz;
            frac_dop = [-.003] * prf_hz;
            v = frac_dop * prf_hz *-wvl/2; % Target velocity (m/s) (causing range walk)

            % Time and range calculations
            num_pulses = 256;
            t = (0:num_pulses-1)/prf_hz; % fast time

            % Calculate the number of samples per pulse
            n_samp_per_pulse  = round(fs_hz/prf_hz);

            % Create the range frequency vector(fast time frequency)
            range_freq = ((0:n_samp_per_pulse-1) * fs_hz/n_samp_per_pulse).';

            % Create range migration slow time vector
            range_migration = R0_m + v(:) * t; % Targets range over pulses
            no_range_migration = R0_m + t;

            % Calculate the phase shift due to range migration i.e.
            % X(-wto);
            phase_shift       = 2 * pi * (fc_hz + range_freq) * (2 * range_migration/c_mps);
            phase_shift_no_rm = 2 * pi * (fc_hz + range_freq) * (2 * no_range_migration/c_mps);

            % Simulate target response
            tx_pulse = rectpuls(t - pulse_width_s/2, pulse_width_s); % Rectangular pulse

            % Transform data into frequency domain
            tx_pulse_freq = fft(tx_pulse,n_samp_per_pulse).';

            % Create fast time( range samples) vs. slow time(PRI's) matrix
            rx_signal_freq = zeros(n_samp_per_pulse, num_pulses);
            rx_signal_freq_no_rm = zeros(n_samp_per_pulse, num_pulses);

            % Apply time-shift to tx pulse in frequency domain
            % Use 2-D phase argument to apply FT time-shift property to
            % shift the signal in time.

            % FT time-shift property states:
            %   A shift of t0 in the time domain is equivalent to
            %   introducing a phase shift of -wt0 in the frequency domain,
            %   but the amplitude stays the same.
            phase_shift_exp = exp(-1j * phase_shift); % exp(–𝝎𝒕0)
            phase_shift_exp_no_rm = exp(-1j * phase_shift_no_rm); % exp(–𝝎𝒕0)
            rx_signal_freq = tx_pulse_freq .* phase_shift_exp; % Apply phase shift(i.e. X(𝝎)* exp(–𝝎𝒕0))
            rx_signal_freq_no_rm = tx_pulse_freq .* phase_shift_exp_no_rm;% Apply phase shift(i.e. X(𝝎)* exp(–𝝎𝒕0))

            % Convert rx signal to time domain
            tmp_data = ifft(rx_signal_freq,[],1);
            tmp_data_no_rm = ifft(rx_signal_freq_no_rm,[],1);

            % Visualize range migration effects on time domain signal
            figure(1)
            subplot(2,1,1)
            tmp_data_mag = abs(tmp_data);
            tmp_data_power_rm = 20*log10(tmp_data_mag);
            imagesc(tmp_data_power_rm);
            title('Time Domain Signal with Range Migration')
            xlabel('Pulse')
            ylabel('Range Samples')
            colormap('jet')
            colorbar;
            axis xy

            subplot(2,1,2)
            tmp_data_mag_no_rm = abs(tmp_data_no_rm);
            tmp_data_power_no_rm = 20*log10(tmp_data_mag_no_rm);
            imagesc(tmp_data_power_no_rm);
            title('Time Domain Signal without Range Migration')
            xlabel('Pulse')
            ylabel('Range Samples')
            colormap('jet')
            colorbar;
            axis xy

            % Verify that range walk exists in the signal by looking at the
            % phase difference between each pulse. it should be
            % approximately constant.

            A = tmp_data_power_no_rm(:,128);
            [max_val,row_idx] = max(A(:));
            [no_rm_row,no_rm_col] = ind2sub(size(A), row_idx);
            fprintf('No RW - Max value: %d at row %d, column %d\n',max_val, no_rm_row, no_rm_col);

            B = tmp_data_power_rm(:,128);
            [max_val,row_idx] = max(B(:));
            [rm_row,rm_col] = ind2sub(size(B), row_idx);
            fprintf('RW - Max value: %d at row %d, column %d\n',max_val, rm_row, rm_col);

            % Check to see if the point with the highest peak in each pulse
            % is at the same index. If the peaks are at the same index then
            % no range migration occured.
            testCase.verifyNotEqual(no_rm_row,rm_row);

        end

        function testDataGenRangeMigratedLFMSignal(testCase)
            % Test: Add range cell migration to a delta function. 
            % This test verifies the range walk efects are added to a 
            % LFM waveform.
            
            % constants
            c_mps = physconst('LightSpeed'); % Speed of light (m/s)

            % Parameters
            fs_hz  = 60e6;              % Sampling frequency (Hz)
            fc_hz  = 10e9;              % Carrier frequency (Hz)
            prf_hz = 50e3;              % Pulse repetition frequency (Hz)
            B_hz   = 10e6;
            pulse_width_s = 1/prf_hz;   % Pulse width (s)
            R0_m = 50;                % Initial range to the target (m)
            lambda_m = c_mps/fc_hz;
            frac_dop = [0.00002] * (prf_hz);
            v = frac_dop * prf_hz * -lambda_m / 2; % Target velocity (m/s) (causing range walk)
            v_max = (prf_hz * -lambda_m)/2;
            fprintf("Max Unambig Velocity: %2.2f mps \n", v_max);
            fprintf("Velocity: %2.2f mps \n", v);

            % Time and range calculations
            n_samples_per_pulse = round(fs_hz/prf_hz);
            n_pulses = 256;        % Number of pulses
            t = (0:n_pulses-1)/prf_hz;   % Fast time (within a pulse)

            % Create the range frequency vector(fast time frequency)
            range_freq = ((0:n_samples_per_pulse-1) * fs_hz/n_samples_per_pulse).';

            % create range walk slow time vector
            range_migration = R0_m + v * t;
            no_range_migration = R0_m + t;

            % Create transmit waveform
            tx_pulse = create_lfm_pulse_samples("Causal",1,"ChirpUpDown",1,"Fs",fs_hz,"NumberOfSamples",n_samples_per_pulse,"SweepBandwidth",B_hz);

            % Transform tx pulse to freq domain and convert to col vector
            tx_pulse_freq  = fft(tx_pulse,n_samples_per_pulse).';
            
            %% Simulate receive pulse with range migration effects

            % Frequency domain phase shift can be used to induce a time
            % shift in the time domain. Also known as the Fourier Transform 
            % shift property.

            % create rx signals
            rx_signal_freq_rm = zeros(n_samples_per_pulse, n_pulses);
            rx_signal_freq_no_rm = zeros(n_samples_per_pulse, n_pulses);

            % Apply time shift using Fourier Transform shift property
            % which states that a shift of t0 in the time domain is
            % equivalent to introducing a phase shift of -wt0 in the
            % frequency domain but the amplitude of the signal remains the
            % same. 

            % Calculate phase shift due to range migration
            phase_shift_no_rm  = 4 * pi * (fc_hz + range_freq) * no_range_migration/c_mps;
            phase_shift_rm     = 4 * pi * (fc_hz + range_freq) * range_migration/c_mps;

            % Apply phase shift
            phase_shift_exp_no_rm = exp(-1j * phase_shift_no_rm);
            rx_signal_freq_no_rm = tx_pulse_freq .* phase_shift_exp_no_rm;

            phase_shift_exp_rm = exp(-1j * phase_shift_rm);
            rx_signal_freq_rm = tx_pulse_freq .* phase_shift_exp_rm;

            % Convert rx signal to fast time vs. slow time to range doppler
            tmp_data_no_rm = ifft(rx_signal_freq_no_rm,[],1);
            tmp_data_no_rm_mf = matched_filter("FftDim", 1,  ...
                 "Signal", tmp_data_no_rm, "Template", tx_pulse_freq);

            if(0)
                % Plot the output signal of the matched filter
                % Plot magnitude response of the matched filter
                figure;
                % Frequency axis
                f = linspace(-fs_hz/2, fs_hz/2, length(tx_pulse));
                plot(f, fftshift(abs(fft(tmp_data_no_rm_mf(:,1)))));
                xlabel('Frequency (Hz)');
                ylabel('Magnitude');
                title('Magnitude Response of Matched Filter');
                grid on;
            end

            % Compute match filter in range dimension
            tmp_data_rm = ifft(rx_signal_freq_rm,[],1);
            tmp_data_rm_mf = matched_filter("FftDim", 1,  ...
               "Signal", tmp_data_rm, "Template", tx_pulse_freq);

            
            %% Now we can visualize range walk effects
            figure(400)
            t = linspace(0, pulse_width_s, fs_hz*pulse_width_s);
            lfm_sig = tx_pulse;
            % extract real an imag parts for analysis
            real_part = real(lfm_sig);
            imag_part = imag(lfm_sig);

            subplot(2,1,1)
            plot(t, real_part,'b', 'LineWidth', 1.5); hold on;
            plot(t, imag_part,'r', 'LineWidth', 1.5); hold off;
            xlabel('Times(s)')
            ylabel('Amplitude');
            title("Real and Imaginary Parts of LFM")
            legend('Real Part','Imaginary Part');
            grid on

            subplot(2,1,2)
            real_part_shifted = real(tmp_data_rm(:,1));
            imag_part_shifted = imag(tmp_data_rm(:,1));
            plot(t, real_part_shifted,'b', 'LineWidth', 1.5); hold on;
            plot(t, imag_part_shifted, 'r', 'LineWidth', 1.5); hold off;
            xlabel('Times(s)')
            ylabel('Amplitude');
            title("Real and Imaginary Parts of Shifted LFM")
            legend('Real Part','Imaginary Part');
            grid on

            figure(401)
            tmp_data_mag_no_rm = abs((tmp_data_no_rm_mf));
            tmp_data_power_no_rm = 20*log10(tmp_data_mag_no_rm);
            imagesc(tmp_data_power_no_rm);
            xlabel('Range Samples')
            ylabel('Pulses');
            axis xy
            colorbar
            title("Fast Time vs. Slow Time no range migration")


            figure(402)
            tmp_data_mag_rm = abs(tmp_data_rm_mf);
            tmp_data_power_rm = 20*log10(tmp_data_mag_rm);
            imagesc(tmp_data_power_rm);
            xlabel('Range Samples')
            ylabel('Pulses');
            axis xy
            colorbar
            title("Fast Time vs. Slow Time w/ range migration")

            % Verify that range walk exists in the signal by looking at the
            % phase difference between each pulse. it should be
            % approximately constant.

            A = tmp_data_power_no_rm(:,128);
            [max_val, row_idx] = max(A(:));
            [no_rm_row, no_rm_col] = ind2sub(size(A), row_idx);
            fprintf('No RW - Max value: %d at row %d, column %d\n',max_val, no_rm_row, no_rm_col);

            B = tmp_data_power_rm(:,128);
            [max_val, row_idx] = max(B(:));
            [rm_row, rm_col] = ind2sub(size(B), row_idx);
            fprintf('RW - Max value: %d at row %d, column %d\n',max_val, rm_row, rm_col);

            % Check to see if the point with the highest peak in each pulse
            % is at the same index. If the peaks are at the same index then
            % no range migration occured.
            testCase.verifyNotEqual(no_rm_row,rm_row);

        end

        function testRangeCellMigrationCorrection(testCase)
            % This function tests that the range cell migration function
            % re-aligns a signal that contains range migration effects
            % constants
            c_mps = physconst('LightSpeed'); % Speed of light (m/s)

            % Parameters
            fs_hz  = 40e6;              % Sampling frequency (Hz)
            fc_hz  = 16e9;              % Carrier frequency (Hz)
            prf_hz = 60e3;              % Pulse repetition frequency (Hz)
            pulse_width_s = 1/prf_hz;
            B_hz   = 20e6;
            R0_m = 1000;                % Initial range to the target (m)
            lambda_m = c_mps/fc_hz;
            frac_dop = [0.00005] * (prf_hz);
            v_mps = frac_dop * prf_hz * -lambda_m / 2; % Target velocity (m/s) (causing range walk)
            v_max = (prf_hz * -lambda_m)/2;
            fprintf("Max Unambig Velocity: %2.2f mps \n", v_max);
            fprintf("Velocity: %2.2f mps \n", v_mps);

            % Time and range calculations
            n_samples_per_pulse = round(fs_hz/prf_hz);
            n_pulses = 2000;        % Number of pulses
            t_pri_s = (0:n_pulses-1)/prf_hz;   % Fast time (within a pulse)

            % Create the range frequency vector(fast time frequency)
            range_freq = ((0:n_samples_per_pulse-1) * fs_hz/n_samples_per_pulse).';

            % Create transmit waveform
            tx_pulse = create_lfm_pulse_samples("Causal",1,"ChirpUpDown",1, ...
                "Fs",fs_hz,"NumberOfSamples",n_samples_per_pulse,"SweepBandwidth",B_hz);
            
            % Transform tx pulse to freq domain and convert to col vector
            tx_pulse_freq  = fft(tx_pulse,n_samples_per_pulse).';
            %% Step 1: Generate signal without migration
            % compute time delay as a function of range
            time_delay = 2 * R0_m / c_mps;
            
            % Compute delay in samples
            sample_delay = time_delay * fs_hz;

            v_const_mps = 1;
            no_range_migration_m = R0_m + v_const_mps * t_pri_s; % meters

            no_range_migration_samples = no_range_migration_m * 2 * fs_hz /c_mps; % samples

            % Calculate phase shift due to range migration
            t0 = 2 * no_range_migration_m / c_mps; % seconds
            phase_shift_no_rm  = 2 * pi * (fc_hz + range_freq) *  t0;
            no_rm_shift_in_samples = t0 * fs_hz;

            % Apply phase shift
            phase_shift_exp_no_rm = exp(-1j * phase_shift_no_rm);
            rx_signal_freq_no_rm = tx_pulse_freq .* phase_shift_exp_no_rm;

            % Convert rx signal to fast time vs. slow time to range doppler
            tmp_data_no_rm = ifft(rx_signal_freq_no_rm,[],1);
            tmp_data_no_rm_mf = matched_filter("FftDim", 1,  ...
                 "Signal", tmp_data_no_rm, "Template", tx_pulse);


            %% Step 2: Generate signal with range migration
            % compute time delay as a function of range
       
            range_migration_m = R0_m + v_mps * t_pri_s;

            range_migration_samples = (range_migration_m * 2 * fs_hz) /c_mps;

            rcm_delta_samples = no_range_migration_samples - range_migration_samples;

            % create rx signals
            rx_signal_freq_rm = zeros(n_samples_per_pulse, n_pulses);

            % Calculate phase shift due to range migration
            t0 = 2 * range_migration_m / c_mps; % seconds
            phase_shift_rm     = 2 * pi * (fc_hz + range_freq) * t0;
            rm_shift_in_samples = t0 * fs_hz;

            rm_shift_in_samples_diff = no_rm_shift_in_samples - rm_shift_in_samples;

            % Apply phase shift
            phase_shift_exp_rm = exp(-1j * phase_shift_rm);
            rx_signal_freq_rm = tx_pulse_freq .* phase_shift_exp_rm;

            % Convert rx signal to fast time vs. slow time to range doppler
            % Compute match filter in range dimension
            tmp_data_rm = ifft(rx_signal_freq_rm,[],1);
            tmp_data_rm_mf = matched_filter("FftDim", 1,  ...
               "Signal", tmp_data_rm, "Template", tx_pulse);

            %% Step 3: Apply range migration correction
            % total duration of pulse(secs)
            % pulse_duration_s = n_samples_per_pulse/fs_hz;
            % % our starting point would be zero instead of number samples/2
            % dt = 1 /fs_hz;
            % time_vec_s = 0:dt:pulse_duration_s-(1/fs_hz);

            % n_samples_per_pulse = round(n_samples_per_pulse  + n_samples_per_pulse/2);
            % temp_rcmc_mat = zeros(n_samples_per_pulse,n_pulses);
            % N = -n_samples_per_pulse/2:n_samples_per_pulse/2-1;
            % 
            % % time vector is number of samples / sampling rate = secs/pulse
            % time_vec_s = N/fs_hz;
            % % for idx = 1:n_pulses
            % %     temp_rcmc_mat(:, idx) = sinc_interp2(tmp_data_rm_mf(:,idx), fs_hz, time_vec_s);
            % % end
            % corrected_data = temp_rcmc_mat;


            %% Plots
            figure(400)
            t = 0:1/fs_hz:pulse_width_s;
            lfm_signal = tx_pulse;
            real_part = real(lfm_signal);
            imag_part = imag(lfm_signal);
            subplot(2,1,1)
            plot(t, real_part, 'b','LineWidth',1.5); hold on;
            plot(t, imag_part, 'r','LineWidth',1.5); hold on;
            xlabel('Time(s)')
            ylabel('Amplitude');
            axis xy
            title("Real and Imaginary Parts of LFM Signal")
            legend('Real Part', 'Imaginary Part');
            grid on

            subplot(2,1,2)
            real_part_shifted = real(tmp_data_rm(:,1));
            imag_part_shifted = imag(tmp_data_rm(:,1));
            plot(t, real_part_shifted, 'b','LineWidth',1.5); hold on;
            plot(t, imag_part_shifted, 'r','LineWidth',1.5); hold on;
            xlabel('Time(s)')
            ylabel('Amplitude');
            axis xy
            title("Real and Imaginary Parts of LFM Signal")
            legend('Real Part', 'Imaginary Part');
            grid on

            figure(401)
            subplot(2,1,1)
            tmp_data_mag_no_rm = abs((tmp_data_no_rm_mf));
            tmp_data_power_no_rm = 20*log10(tmp_data_mag_no_rm);
            imagesc(tmp_data_power_no_rm);
            xlabel('Range Samples')
            ylabel('Pulses');
            axis xy
            colorbar
            title("Fast Time vs. Slow Time no range migration")

            subplot(2,1,2)
            tmp_data_mag_rm = abs(((tmp_data_rm_mf)));
            tmp_data_power_rm = 20*log10(tmp_data_mag_rm);
            imagesc(tmp_data_power_rm);
            xlabel('Range Samples')
            ylabel('Pulses');
            axis xy
            colorbar
            title("Fast Time vs. Slow Time w/ range migration")

            % figure(403)
            % tmp_data_mag_rm = abs(corrected_data);
            % tmp_data_power_rm = 20*log10(tmp_data_mag_rm);
            % imagesc(tmp_data_power_rm.');
            % xlabel('Range Samples')
            % ylabel('Pulses');
            % axis xy
            % colorbar
            % title("Fast Time vs. Slow Time w/ RCMC")

            %% Step 4: Verify results
            % Measure residual migration
            % residual_migration = RangeCellMigrationCorrectionTest.estimate_range_migration(corrected_data, no_rm_shift_in_samples, rm_shift_in_samples_diff);
            % 
            % % Define an acceptable error tolerance
            % tolerance = 0.5;
            % 
            % % Verify correction
            % testCase.verifyLessThan(max(abs(residual_migration)), tolerance, ...
            %     "RCMC failed: Residual migration exceeds tolerance.");
        end
    end
end