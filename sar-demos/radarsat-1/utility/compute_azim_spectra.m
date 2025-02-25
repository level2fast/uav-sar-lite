% compute_azim_spectra.m
% ----------------------
%
% Divides the range swath (range cells) of this block into Nspect sections 
% and finds the average power spectrum for each segment.
%
% Fits a parabola to the baseband Doppler centroid vs. range
% This fit is most useful if the whole range swath is selected, but the
% computing times get large if more than 1024 range lines (azimuth cells)
% are selected in one run.
%
% Run "specify_run_parameters.m" and "extract_data.m" first to extract
% the data from the CD. The data can be stored in MAT or binary files.
% -------------------------------------------------------------------------

%  Load the parameters for this run
clear,    format compact
plot_az_specta_results = 0;
set( 0, 'DefaultTextFontSize',   12 )  % Plotting defaults
set( 0, 'DefaultLineLineWidth', 1.5 )
set( 0, 'DefaultAxesFontSize',    8 )
tfs = 13;    lfs = 11;
load CD_run_params

disp ' '
disp '---------------------------------------------------------'
fprintf(' UBC RRSG - Plot the azimuth spectrum of each data block')
disp ' '
disp '---------------------------------------------------------'

Nrowsg = 3;     % Number of subplots in row direction of the figure 
Ncolsg = 3;     % Number of subplots in column direction of the figure 
Nspect = Nrowsg*Ncolsg;              % Total number of spectra to calculate
Nrglpb = floor( Nrg_cells/Nspect );  % No. of range cells per spectra
wd = 0.81/Ncolsg;   dx = wd + 0.045;   x0 = 0.07;    
ht = 0.39/Nrowsg;   dy = 0.28;         y0 = 1-dy;  

for b = 1 : Nblocks % loop through each range block,  2 total are used
    
    file_pre = strcat( output_path, output_prefix, '_', num2str(b) );
    
    disp ' '
    disp (['Load or Extract AGC setting and Data for block ' num2str(b) ])

    % Load a block of 'AGC_values'
    % These values are used to maintain a consistent signal level despite
    % varying return signal strengths caused by different target distances,
    % sizes, reflectivity. It's  primary purpose is improving SNR. AGC helps
    % to maintain an optimal SNR by adjusting the gain to suit the current 
    % signal conditions.
    AGC_values = load_AGC_block( file_pre, first_rg_line, ...
                                          Nrg_lines_blk, b , UseMATfiles );
                                      
    %  Load a block of raw SAR data
    data = load_DATA_block( file_pre, output_path, Nrg_lines_blk, ...
                             Nrg_cells, AGC_values, b, UseMATfiles );
        
    disp 'Compute and plot azimuth power spectra'
    tic
    if(plot_az_specta_results)
        figure, clf
    end

    % Frequency vector for each range line divided into PRFs per rng line
    freq = [0:Nrg_lines_blk-1]*PRF/Nrg_lines_blk;

    % Beginning of the first range line for current block
    b1 = first_rg_line + (b-1) * Nrg_lines_blk;

    % End of the last range line for current block
    b2 = first_rg_line + b * Nrg_lines_blk - 1;
    
    %  Loop through each spectra
    for range_line_idx = 1 : Nspect 

        % Calculate indexes of range cells that will be used to find the az
        % power spectra to be averaged( r1 = range cell 1, r2 = range cell 2
        % there are Nrglpb range cells per spectra
        r1 = 1 + (range_line_idx-1)*Nrglpb;   r2 = r1 + Nrglpb - 1;

        % Calculate the fft of column of data from start r1 to 
        % r2 which is a total of 227 range cells.
        DATA = fft( data(:,r1:r2) );
        
        %TODO SDD: save original range lines used for calculating the average
        % power of the azimuth spectra for testing. This data can be used 
        % to verify that the baseband dopp function produces the correct 
        % output. The rows represent the data extracted from the CD at the 
        % 1st range line. The columns represent the data extracted from the
        % CD starting with the 1st range cell
        if(b==1)
            test_data_for_est_baseband_dopp(:,:,range_line_idx) = data(:,r1:r2);
        end

        % Calculate the average power across range cell 1 to range cell 2
        DATA_TRANSP = DATA.';
        % temp2 = mean(temp1);
        % temp3 = temp2/10e5;
        DATA_aver(range_line_idx,:) = mean( abs(DATA_TRANSP).^2 )/1000000;

        % Scale the average power values from range cell 1 to range cell 2
        ysc(range_line_idx) = 1.05*max(DATA_aver(range_line_idx,:));

    end  % of for krg = 1 : Nspect
    ysc0 = max( ysc );      %  Common vertical scaling for all the spectra
    
    % Loop through the azimuth spectra
    for range_line_idx = 1 : Nspect 

        if(plot_az_specta_results)
            subplot(Nrowsg, Ncolsg, range_line_idx)
            % Plot the average power spectra with freq as x axis
            plot( freq, DATA_aver(range_line_idx,:) ),   grid,   hold on
            set( gca, 'Pos',...
               [x0+dx*mod((range_line_idx-1),Ncolsg)  y0-dy*floor((range_line_idx-1)/Ncolsg)  wd ht])
            axis([0 PRF  0 ysc0]);
        end
        
        % Calculate fft on each column of averaged spectra, so each
        % averaged sub-block to get the azimuth/doppler spectrum
        azim_spec = fft( DATA_aver(range_line_idx,:) )/ Nrg_lines_blk;

        % From Pg. 506 of digital processing of synthetic aperture radar
        % data. "When a sine wave model is used, the estimation of the
        % spectral peak can be found simply as the phase angle of the
        % fundamental harmonic of the magnitude spectrum(i.e. the phase
        % angle of the DFT coefficient corresponding to the one cycle per
        % record)."
        % Get phase angle in radians of 2nd index in the Average power azimuth 
        % spectra in matrix. This corresponds to the fundamental harmonic, 
        % which corresponds to the first non-zero frequency component after the 
        % DC component (which is at index 1 in MATLAB, corresponding to 0 Hz). 
        % The 1st(fundamental) harmonic angle will be negative so multipy by
        % (-1). The angle of first harmonic describes the amount by which 
        % this sample is out of phase with the fundamental frequency. That
        % value corresponds to an estimate of the percentage of the PRF for
        % this sub-block, which is the baseband doppler frequency(i.e the 
        % estimation of the freuqency of the spectral peak. pg. 506)
        angle_first_harmonic = -angle( azim_spec(2) );

        % Calcuate Fdc(Baseband Doppler Centroid) for the current row using
        % the phase angle. 
        perc =angle_first_harmonic / (2*pi);
        Ffrac(range_line_idx) = PRF * angle_first_harmonic / (2*pi);

        baseband_dopp_cent_solution(b,range_line_idx) = Ffrac(range_line_idx);

        % if baseband doppler frequency is in negative part of spectrum add
        % a full PRF make to get the positive value. This works because
        % doppler frequency is alwasy between -PRF/2 to PRF/2.
        if Ffrac(range_line_idx) < 0   
            Ffrac(range_line_idx) = Ffrac(range_line_idx) + PRF;   
        end

        % Fit azimuth spectrum to a sine wave
        sine_fit = real(azim_spec(2)) * cos(2*pi*freq/PRF) - ...
                imag(azim_spec(2)) * sin(2*pi*freq/PRF) + 0.5*azim_spec(1);
        if(plot_az_specta_results)
            plot( freq, 2*sine_fit, 'r--' ); 
        end
       
        % Select start and end range cells for current spectra and block
        r1 = 1 + (range_line_idx-1)*Nrglpb;   r2 = r1 + Nrglpb - 1;



        % plot aizmuth spectra
        if(plot_az_specta_results)
            title(sprintf('RC %4d - %4d   Fdc =%6.0f',...
                r1+first_rg_cell-1, r2+first_rg_cell-1, Ffrac(range_line_idx) ), ...
                'FontSize', lfs );
            if range_line_idx > Nspect - Ncolsg
                xlabel('Azimuth frequency  (Hz)  \rightarrow', 'FontSize', lfs )
            end
            if mod(range_line_idx,Ncolsg) == 1
                ylabel('Power  \rightarrow', 'FontSize', lfs )
            end
            if range_line_idx == 1
                text( 1.55*PRF, 1.7*double(ysc0), sprintf(...
                 'Azimuth power spectra of range lines%6.0f  to%6.0f  and sine fit',...
                  first_rg_line+(b-1)*Nrg_lines_blk, first_rg_line+b*Nrg_lines_blk-1 ),...
                  'Hor', 'center', 'FontSize', tfs )
            end
        end
        pause(0.1)
    end  % of for krg = 1 : Nspect
    toc
    
    %  Plot Ffrac vs. range
    if(plot_az_specta_results)
        figure(203),   clf
        plot( Ffrac, 'bx-', 'MarkerS', 9 ),   grid,   hold on
        eval(['save Ffrac_' num2str(first_rg_line) ' Ffrac'])
        axis([0.5 Nspect+0.5  200 800])
        xlabel('Range swath number', 'FontSize', lfs+1 )
        ylabel('Frequency  (Hz)', 'FontSize', lfs+1 )
        title(sprintf(...
            'Baseband Doppler centroid over%5.0f lines starting at%6.0f',...
            Nrg_lines_blk, first_rg_line+(b-1)*Nrg_lines_blk ), 'FontSize', tfs )
    end
    
    % Find coefficients of Baseband Doppler centroid 2nd order polynomial
    coeff = polyfit( [1:Nspect], Ffrac, 2 );
    if(b == 1)
        % Evaluate the coefficient polynomial for each spectra 1-9
        FitBlk1 = polyval( coeff, [1:Nspect] );
        Fit = FitBlk1;
    elseif(b==2)
        % Evaluate the coefficient polynomial for each spectra 1-9
        FitBlk2 = polyval( coeff, [1:Nspect] );
        Fit = FitBlk2;
    end
    
%     % Evaluate the coefficient polynomial for each spectra 1-9
%     Fit = polyval( coeff, [1:Nspect] );
    % Plot Frequency fit for for each block of range cells
    if(plot_az_specta_results)
        plot( Fit, 'r-' )
        text( 0.79*Nspect, 750, sprintf('C_2 =%6.2f', coeff(1) ), 'FontSize',13 )
        text( 0.79*Nspect, 650, sprintf('C_1 =%7.1f', coeff(2) ), 'FontSize',13 )
        text( 0.79*Nspect, 550, sprintf('C_0 =%5.0f', coeff(3) ), 'FontSize',13 )
        text( 0.13*Nspect, 250, sprintf(...
            'Range cells %3.0f to %4.0f', first_rg_cell, ...
             first_rg_cell+Nrg_cells-1 ), 'FontSize',13 )
    end

    pause(0.1)
  % Save the plot as .eps file 
  % file_eps=strcat(output_path,output_prefix,'azimuth_',num2str(b),'.eps')
  % saveas(gcf,file_eps,eps)
end  % of for b = 1 : Nblocks
beep,   pause(0.3),   beep,   pause(0.3),   beep