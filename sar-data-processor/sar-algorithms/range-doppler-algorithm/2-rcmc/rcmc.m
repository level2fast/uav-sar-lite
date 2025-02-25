function [xr] = rcmc(data, Fn, lambda, vEff, R0, Fs, time_axis_s)
%RCMC
% Range Cell Migration Correction, which is range time and azimuth frequency
% dependent, is performed in the range Doppler domain, where a family of
% trajectories at the same range are transformed into one single
% trajectory. RCMC straightens out these trajectories so that they now run
% parallel to the azimuth frequency axis.
% 
% Inputs:
% 
% Outputs:
%

% Intialize deltaR(fn). This is the amount of range cell migration to
% correct for and represents the target displacement as a function of
% azimuth frequeny fn.

xn      = data; % sampled signal
nIdx    = 1;    % sample counter
t_s     = time_axis_s;   % new time axis sampled at higher Fs
window  = 1;    % window function to cutoff sinc

sincTrain = zeros(length(t_s), length(xn)); % sinc train for each sample

for n = -floor(length(xn)/2):floor(length(xn)/2) % for each sample in the sampled signal
    % Step 1: Calculate the amount range migration in time
    deltaR_fn = calc_rcm_samples(lambda, R0, Fn, vEff, Fs);
        
    %step 2: Calculate a sinc train shifted by amount RCM
    nShift = n-deltaR_fn; % sample shift for this sample
    sincKernel = sinc((t_s - nShift*Ts)/Ts);

    %step 3: Calculate windowed sinc 
    windowed_sinc = sincKernel * window;

    %step 4: Scale sinc windowed sinc by sample x[n] and store result
    sincTrain(:, nIdx) = xn(nIdx) * windowed_sinc;

    % increment counter
    nIdx = nIdx + 1;
end

xr = sum(sincTrain, 2);
end

