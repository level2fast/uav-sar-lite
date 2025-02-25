classdef RngDoppAlgTests < matlab.unittest.TestCase
    properties
        plot_sample_data logical = false;
        slowTime_s
        numpulses
        truncrangesamples
        radarPlatform
        transmitter
        radiator
        collector
        receiver
        channel
        target
        pointTargets
        targetpos_m
        targetvel_mps
        rxsig
        waveform
        flightDuration_s
        tpd_s
        fs_hz
        fastTime
        prf_hz
    end
    methods(TestClassSetup)
        % Shared setup for the entire test class
        function setup_sim_data(testClass)
            %% generate full path
            % Get the full path to the current script
            currentFile = mfilename('fullpath');
            
            % Get the directory of the current script
            [currentFolder, ~, ~] = fileparts(currentFile);
            
            % Define the relative path to the target file from the current script's directory
            relativePath = fullfile(currentFolder, '../', 'common');  % Adjust as needed
            addpath(relativePath);
            
            %% setup test scenario
            c_mps = physconst('LightSpeed');
            fc_hz = 4e9;
            rangeResolution = 3;  
            crossRangeResolution = 3;
            bw_hz = c_mps/(2*rangeResolution);
            testClass.prf_hz = 1000; 
            aperture = 4;  
            testClass.tpd_s = 3*10^-6; 
            testClass.fs_hz = 120*10^6;
            testClass.waveform = phased.LinearFMWaveform( ...
            'SampleRate', testClass.fs_hz, ...
            'PulseWidth', testClass.tpd_s, ...
            'PRF', testClass.prf_hz,...
            'SweepBandwidth', bw_hz);
            speed_mps = 100;  
            testClass.flightDuration_s = 4;
            testClass.radarPlatform  = phased.Platform('InitialPosition', [0;-200;500], 'Velocity', [0; speed_mps; 0]);
            testClass.slowTime_s = 1/testClass.prf_hz;
            testClass.numpulses = testClass.flightDuration_s/testClass.slowTime_s + 1;
            
            maxRange_m = 2500;
            testClass.truncrangesamples = ceil((2*maxRange_m/c_mps)*testClass.fs_hz);
            testClass.fastTime = (0:1/testClass.fs_hz:(testClass.truncrangesamples-1)/testClass.fs_hz);
            % Set the reference range for the cross-range processing.
            antenna = phased.CosineAntennaElement('FrequencyRange', [1e9 6e9]);
            antennaGain = aperture2gain(aperture,c_mps/fc_hz); 
            
            testClass.transmitter = phased.Transmitter('PeakPower', 50e3, 'Gain', antennaGain);
            testClass.radiator = phased.Radiator('Sensor', antenna,'OperatingFrequency', fc_hz, 'PropagationSpeed', c_mps);
            
            testClass.collector = phased.Collector('Sensor', antenna, 'PropagationSpeed', c_mps,'OperatingFrequency', fc_hz);
            testClass.receiver = phased.ReceiverPreamp('SampleRate', testClass.fs_hz, 'NoiseFigure', 30);

            testClass.channel = phased.FreeSpace('PropagationSpeed', c_mps, 'OperatingFrequency', fc_hz,'SampleRate', testClass.fs_hz,...
                'TwoWayPropagation', true);

            testClass.targetpos_m = [800,0,0;1000,0,0; 1300,0,0]'; 
            %testClass.targetpos_m = [800,0,0]'; 
            
            testClass.targetvel_mps = [0,0,0;0,0,0; 0,0,0]';
            %testClass.targetvel_mps = [0,0,0]';
            
            testClass.target = phased.RadarTarget('OperatingFrequency', fc_hz, 'MeanRCS', [10,10,10]);
            %testClass.target = phased.RadarTarget('OperatingFrequency', fc_hz, 'MeanRCS', [10]);

            testClass.pointTargets = phased.Platform('InitialPosition', testClass.targetpos_m,'Velocity',testClass.targetvel_mps);
            if(testClass.plot_sample_data)
                % The figure below describes the ground truth based on the target
                % locations.
                figure;h = axes;
                plot(testClass.targetpos_m(2,1),testClass.targetpos_m(1,1),'*g');
                hold all;
                plot(testClass.targetpos_m(2,2),testClass.targetpos_m(1,2),'*r');
                hold all;
                plot(testClass.targetpos_m(2,3),testClass.targetpos_m(1,3),'*b');
                hold off;
                set(h,'Ydir','reverse');
                xlim([-10 10]);
                ylim([700 1500]);
                title('Ground Truth');
                ylabel('Range');
                xlabel('Cross-Range');
            end

        end
        function perform_data_collection(testClass)
            % Define the broadside angle
            refangle = zeros(1,size(testClass.targetpos_m,2));
            testClass.rxsig = zeros(testClass.truncrangesamples,testClass.numpulses);
            for idx = 1:testClass.numpulses
                % Update radar platform and target position
                [radarpos_m, radarvel] = testClass.radarPlatform(testClass.slowTime_s);
                [targetpos, targetvel] = testClass.pointTargets(testClass.slowTime_s);
                
                % Get the range and angle to the point targets
                [targetRange_m, targetAngle_rad] = rangeangle(targetpos, radarpos_m);
                
                % Generate the LFM pulse
                sig = testClass.waveform();

                % Use only the pulse length that will cover the targets.
                sig = sig(1:testClass.truncrangesamples);
                
                % Transmit the pulse
                sig = testClass.transmitter(sig);
                
                % Define no tilting of beam in azimuth direction
                targetAngle_rad(1,:) = refangle;
                
                % Radiate the pulse towards the targets
                sig = testClass.radiator(sig, targetAngle_rad);
                
                % Propagate the pulse to the point targets in free space
                sig = testClass.channel(sig, radarpos_m, targetpos, radarvel, targetvel);
                
                % Reflect the pulse off the targets
                sig = testClass.target(sig);
                
                % Collect the reflected pulses at the antenna
                sig = testClass.collector(sig, targetAngle_rad);
                
                % Receive the signal  
                testClass.rxsig(:,idx) = testClass.receiver(sig);
                
            end
        end
    end

    methods(Test)
        % Test methods
        function calcRngOfClosestApproachTest(testCase)
            % This function should find the smallest distance among all possible
            % ranges between the platform and each target when given a test 
            % data set with targets at varying ranges. The range of closest
            % approach acts as a reference point to help align phase
            % measurements for SAR in processing steps further down the
            % pipeline. 
            
            plot_data = false;
            % To calculate th range of closest approach I need a sample set of
            % SAR data and the following params: slant range,
            % velocity of radar, and slow time vector. 

            signal_data = testCase.rxsig;
            %signal_data = (pulsint(testCase.rxsig,'noncoherent'));
            fft_size = 2^(nextpow2(size(signal_data,1)));
            tx_waveform_samples = testCase.waveform();

            % calculate the slant range to all targets
            slant_range_vec = matched_filter(Signal=signal_data, ...
                Template=tx_waveform_samples(1:testCase.truncrangesamples), ...
                LenFFT=fft_size, ...
                FftDim=1, ...
                plot_signal=true);

            % Find the number of range bins and azimuth positions
            [num_range_bins, num_azimuth_positions] = size(slant_range_vec(1:testCase.truncrangesamples,:));

            %% find the target with the greatest magnitude for this pulse
            max_peak_idx = calc_max_peak_from_data(slant_range_vec, num_range_bins);

            % Calculate range of closest appraching using the peak with the 
            % highest amplitude in the compressed signal.
            wfrm = testCase.waveform;
            bw = wfrm.SweepBandwidth;
            %range_axis = (0:num_range_bins - 1) * (physconst('LightSpeed') / (2 * bw));
            %max_peak_rng = range_axis(max_peak_idx);  % Corresponding range
            
            velocity_mps = max(testCase.radarPlatform.Velocity);
            slow_time_s = (0:testCase.numpulses - 1) /testCase.prf_hz;
 
            rng_of_closest_appr = calc_rng_of_closest_approach( ...
                slant_range_m=max_peak_idx, ...
                slow_time_secs=slow_time_s, ...
                vr_mps=velocity_mps);

            expected = 800; % meters
            testCase.verifyEqual(rng_of_closest_appr, expected, 'AbsTol', 70)

            if(plot_data)
                figure
                imagesc(real(testCase.rxsig));
                title('SAR Raw Data')
                xlabel('Cross-Range Samples')
                ylabel('Range Samples')

                figure
                imagesc(real(slant_range_vec));
                title('SAR Range Compressed Data')
                xlabel('Cross-Range Samples')
                ylabel('Range Samples')
            end

        end
       
    end
end