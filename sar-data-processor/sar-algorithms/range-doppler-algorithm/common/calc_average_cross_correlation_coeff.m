function [accc_val] = calc_average_cross_correlation_coeff(AcccParams)
%CALC_AVERAGE_CROSS_CORRELATION_COEFF Calculate the average cross
%correlation coefficient for 1 range line of data.
% 
% The conjugate product between samples must be calculated to find the 
% angle of the vector that represents the phase difference between each 
% sample. The sum can be taken instead of the average, since only the angle
% is of interest. This estimate is called the "Average Cross Correlation Coefficient(ACCC)"
% at lag one, because the sum is the component of the autocorrelation function of the
% signal when the shift is one sample. This creates a vector whose angle can
% be used to approximate the average phase change
% ACCC = sum(s*(n) .* s(n + Δn)) where s(n) = odd samples and 
% s(n + Δn) = even samples.
%
%   Inputs
%       slow_time_data - range compressed data
%
%   Outputs
%       ACCC - a vector whose angle can be used to approximate the average
%       phase change of 1 range line(vector) of data

arguments
    AcccParams.slow_time_data (1,:) {mustBeVector} = []
end
data = AcccParams.slow_time_data;

% First grab even and odd indices which correspond to s(n) & s(n + 1) for
% each azimuth sample
data_len = length(data);
is_even = mod(data_len,2);
sig_inc_vec_sum = complex(0,0);

% next check to see if the # of samples is even or odd
if(is_even == 0 )
    odd_indices = 1:2:length(data);
    odd_elements = data(odd_indices);
    
    even_indices = 2:2:length(data);
    even_elements = data(even_indices);
    
    % Compute the average cross correlation coefficient by peforming an
    % element wise multiplication.This effectively multiplies each sample
    % in the azimuth signal with the next sample since the 
    % ACCC = s*(n) .* s(n + Δn)
    sig_inc_vec_sum = sum(conj(odd_elements) .* even_elements); 
else % odd # of samples
    for idx = 1:(data_len - 1)
        sig_inc_vec = [conj(data(idx)) * data(idx+1), sig_inc_vec_sum];
        sig_inc_vec_sum = sum(sig_inc_vec,"all");     
    end
end

% Calculate the angle of Average Cross Correlation Coeffeicient
accc_val = sig_inc_vec_sum;
end

