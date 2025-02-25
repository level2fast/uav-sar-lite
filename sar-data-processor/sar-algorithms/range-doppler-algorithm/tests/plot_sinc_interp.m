function [reconstructed_signal] = plot_sinc_interp(plotSincInterp)
%PLOT_SINC_INTERP Summary of this function goes here
%   Detailed explanation goes here
arguments
plotSincInterp.time_axis;
plotSincInterp.cont_signal;
plotSincInterp.sampling_period;
plotSincInterp.sampled_time_axis;
plotSincInterp.sample_signal;
end
t = plotSincInterp.time_axis;
xc = plotSincInterp.cont_signal;
Ts = plotSincInterp.sampling_period;
tn = plotSincInterp.sampled_time_axis;
xn = plotSincInterp.sample_signal;

% Plot the CT signal and sampled signal
figure
hold on
grid on
plot(t, xc)
stem(tn, xn, 'o')
legend('"Continuous time signal"', 'Sampled signal')
xlabel('t(secs)')
ylabel('Amplitude')

% Create and plot sinc train
sincTrain = zeros(length(t), length(xn));
nind = 1;
figure
cmap = colormap(jet(length(-floor(length(xn)/2):floor(length(xn)/2))));
ax = axes('colororder', cmap);
hold on
grid on

plot(t, xc, 'k', 'LineWidth', 3)
for n = -floor(length(xn)/2):floor(length(xn)/2)
    
   sincTrain(:, nind) = xn(nind)*sinc((t - n*Ts)/Ts);
   p = plot(t, sincTrain(:, nind), 'LineWidth', 2);
   stem(tn(nind), xn(nind), 'Color', p.Color, 'LineWidth', 2)
   nind = nind + 1;
end
xlabel('t')
ylabel('Amplitude')
set(gca, 'FontSize', 20, 'LineWidth', 3, 'FontWeight', 'bold')

xr = sum(sincTrain, 2); % sum(sincTrain, 2) is the interpolated/reconstructed signal, should be equal to xc
rmse = mean(abs(xc - xr).^2);
figure
hold on
grid on
plot(t, xc)
plot(t, xr) 
txt = {'RMSE:', num2str(rmse)};
text(0,0,txt);
legend('Original Time Signal', 'Reconstructed Time signal')
title("Continuous & Reconstructed Signals")
hold off
reconstructed_signal = xr;
end

