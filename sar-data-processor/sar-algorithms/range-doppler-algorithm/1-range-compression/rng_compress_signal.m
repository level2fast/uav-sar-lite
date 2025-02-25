function [compressed_signal] = rng_compress_signal(rngCompSig)
% compress_signal: Function to compress an Linear FM chirp signal
%   This function takes a Linear FM Chirp and a signal with noise
%   as input and computes the match filter to be returned to the caller.
%   This function takes the complex conjugate of the DFT of the zero padded
%   pulse replica length(minus one), the FFT length should be a few times longer
%   than the pulse replica length to obtain efficient processing.
%   corresponds to option 2 on page 95 of "Digital Processing of Synthetic
%   Aperture Radar Data" by Ian G. Cumming and Frank H. Wong. ISBN-13: 978-1-58053-058-3.
arguments
    rngCompSig.ref_sig           (1,:) {mustBeVector} = []
    rngCompSig.recv_sig          (:,:)  = []
    rngCompSig.perform_ifft      (1,1) {mustBeInRange(rngCompSig.perform_ifft,0,1)} = 0
    rngCompSig.number_of_samples (1,1) {mustBeInteger} = 0
end

reference_signal  = rngCompSig.ref_sig;
received_signal   = rngCompSig.recv_sig;
peform_ifft       = rngCompSig.perform_ifft;
number_of_samples = rngCompSig.number_of_samples;

[rows, ~] = size(received_signal);
fft_size = (number_of_samples) *2;

% Transform range to Frequency domain by performing an
% FFT on across range cells(i.e. row) of data. See FFT(X,[],DIM) or fft(X,n,2) in
% this case returns the n-point Fourier transform of each row.
noisy_signal_freq_dm = fft((received_signal),fft_size,2);

% Compute fft on reference signal
reference_signal_freq_dm = fft(reference_signal,fft_size,2);

% reference_signal_freq_dm = repmat(reference_signal_freq_dm, [rows 1]);

% Compute Matched filter.
restored_signal = noisy_signal_freq_dm .* conj(reference_signal_freq_dm);

% Perform range IFFT
if(peform_ifft)

    compressed_signal = ifft(restored_signal);
    compressed_signal(:,size(reference_signal,2)+1:end) = []; % throwaway region

else
    compressed_signal = restored_signal;
end

end

