classdef SincInterpolationTest < matlab.unittest.TestCase
    properties
        plot_figures = false
    end

    methods(TestClassSetup)
        % Shared setup for the entire test class
    end
    
    methods(TestMethodSetup)
        % Setup for each test
    end
    
    methods(Test)
        % Test methods
        function testSincInterp2Function(testCase)
            % Create "continuous time" signal, Fc >> f
            FS_cont_hz = 1e6; % very high sample rate to simulate "continuous time"
            Tc_secs = 1/FS_cont_hz; % sampling period
            t_secs = (-0.1:Tc_secs:0.1)'; % time axis
            f_hz = 10; % signal frequency
            continuous_signal = sin(2*pi*f_hz*t_secs); % "continuous time" signal
            
            % Create sampled signal
            Fs_hz = 100; % sampling rate
            Ts_secs = 1/Fs_hz; % sampling period
            ratio = round(Ts_secs/Tc_secs);
            tn_secs = t_secs(1:ratio:end); % sampled time axis
            x_n = continuous_signal(1:ratio:end); % sampled signal
            
            % Create and plot sinc train
            sincTrain = zeros(length(t_secs), length(x_n));
            n_ind = 1;
            if(testCase.plot_figures)
                % Plot the CT signal and sampled signal
                figure
                hold on
                grid on
                plot(t_secs, continuous_signal)
                stem(tn_secs, x_n, 'o')
                legend('"Continuous time signal"', 'Sampled signal')
                xlabel('t(secs)')
                ylabel('Amplitude')

                figure
                cmap = colormap(jet(length(-floor(length(x_n)/2):floor(length(x_n)/2))));
                ax = axes('colororder', cmap);
                hold on
                grid on
                
                plot(t_secs, continuous_signal, 'k', 'LineWidth', 3)

                for n = -floor(length(x_n)/2):floor(length(x_n)/2)
                   sincTrain(:, n_ind) = x_n(n_ind)*sinc((t_secs - n*Ts_secs)/Ts_secs);
                   p = plot(t_secs, sincTrain(:, n_ind), 'LineWidth', 2);
                   stem(tn_secs(n_ind), x_n(n_ind), 'Color', p.Color, 'LineWidth', 2)
                   n_ind = n_ind + 1;
                end

                % sum(sincTrain, 2) is the interpolated/reconstructed signal, should be equal to xc
                reconstructed_signal = sum(sincTrain, 2); 
                rmse = mean(abs(continuous_signal - reconstructed_signal).^2);    
                xlabel('t(secs)')
                ylabel('Amplitude')
                set(gca, 'FontSize', 20, 'LineWidth', 3, 'FontWeight', 'bold')
                figure
                xlabel('t(secs)')
                ylabel('Amplitude')
                hold on
                grid on
                plot(t_secs, continuous_signal)
                plot(t_secs, reconstructed_signal)
                txt = {'RMSE:', num2str(rmse)};
                text(0,0,txt);
                legend('Original Time Signal', 'Reconstructed Time signal')
                title("Continuous & Reconstructed Signals")
                hold off
            else
                T = 0.1;
                dt = 1/FS_cont_hz;
                new_time_axis = -T:dt:T;
                [reconstructed_signal] = sinc_interp2(x_n,Fs_hz,new_time_axis);
            end

            rmse = mean(abs(continuous_signal - reconstructed_signal).^2);
            rmse_treshold = 1e-3;
            testCase.verifyLessThan(rmse, rmse_treshold);
        end
        function testSincInterpFunction(testCase)
            % Sinc interpolation is a foundational tool that can be used to
            % perform range cell migration correction. This function tests
            % sinc interpolation of basic sine wave
            Fs = 10;   % Original sampling frequency
            Ts = 1/Fs; % sampling period of original signal
            t = 0:Ts:1; % Original sampling times
            f = 10;
            orig_sig = sin(2 * pi * f * t); % original signal

            % make interpolated time series
            newFs = 100;
            newTs = 1/newFs;
            t_interp = 0:newTs:1; % New time vector (finer resolution)

            % Expected result
            expected_sig = sin(2 * pi * f * t_interp);

            % Sinc interpolation
            x_interp = sinc_interp("Ts_s",Ts,"interp_time_axis_s",t_interp,"time_axis_s",t,"sampled_signal",orig_sig);

            if(testCase.plot_figures)
                % Plot results
                figure;
                subplot(2, 1, 1);
                stem(t, orig_sig, 'r', 'LineWidth', 1.5); % Original sampled points
                hold on;
                plot(t_interp, x_interp, 'b');      % Interpolated signal
                xlabel('Time (s)');
                ylabel('Amplitude');
                title('Sinc Interpolation of a Sine Wave');
                legend('Original Samples', 'Interpolated Signal');
                grid on;
                
                subplot(2, 1, 2);
                plot(t_interp, sin(2 * pi * f * t_interp), 'k--'); % True continuous signal
                hold on;
                plot(t_interp, expected_sig, 'b');       % Interpolated signal
                xlabel('Time (s)');
                ylabel('Amplitude');
                title('Comparison with True Continuous Signal');
                legend('True Signal', 'Interpolated Signal');
                grid on;
            end

            % calculate RMSE
            rmse = mean(abs(expected_sig - x_interp).^2);
            % Assert interpolation is accurate
            testCase.verifyLessThanOrEqual(rmse, .001);
        end
    end
    
end