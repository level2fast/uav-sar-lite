function output = matched_filter(matchFilter)
    % matched_filter - Performs matched filtering of an input signal with a template
    %
    % Syntax:
    % output = matched_filter(signal, template)
    %
    % Inputs:
    %   signal - The input signal to be filtered (vector)
    %   template - The template (reference signal) for matched filtering (vector)
    %
    % Outputs:
    %   output - The result of the matched filter operation (vector)
    %
    % Example:
    %   output = matched_filter(signal, template);
    %   plot(output);
    arguments
        matchFilter.Signal   (:,:) {mustBeNonempty} = []
        matchFilter.Template (:,:) {mustBeNonempty} = []
        matchFilter.FftDim   (1,1) {mustBeNonempty} = 0
    end
    signal   = matchFilter.Signal;
    template = matchFilter.Template;
    dim      = matchFilter.FftDim;
    len_signal = size(template,2);
    len_fft = size(signal,1) + size(template,1) - 1;

    signal_freq   = fft(signal, len_fft, dim);
    template_freq = fft(template(:), len_fft);

    % Take the complex conjugate of the template (matched filter operation)
    matched_template = conj(template_freq);
       
    output_freq = matched_template.* signal_freq;
    output = ifft(output_freq, [], 1);
    output = output(1:len_signal,:);

end
