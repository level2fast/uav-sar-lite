classdef RangeDopplerAlg < SarAlgorithmIf 
    %RangeDopplerAlg Provides steps for performing Range Doppler Algorithm
    %   This class performs all SAR processing using the Basic Range Doppler
    %   Algorithm. The steps involved in this processing were taken from the book 
    %   "Digital Processing of Synthetic Aperture Radar Data" by Ian G. Cumming and 
    %   Frank H. Wong. ISBN-13: 978-1-58053-058-3. 
    %   The RDA will be performed using the following steps:
    %   Range Compression -> Azimuth FFT -> RCMC -> Azimuth Compression
    %   -> Azimuth IFFT and Look Summation->Compressed Data
    % 
    %   Each step is broken down into sub-steps which contain another set of
    %   steps specific to it. This program is written in a modular way to allow for
    %   easy reuse of certain parts of the code if needed. 
    properties(Access=protected)
        % Data Processor Params
        Operation       {mustBeText,mustBeNonempty,mustBeMember(Operation,{'Low Squint','Exact SRC','Approximate SRC'})} 
        Multilooking    {mustBeText,mustBeNonempty,mustBeUnderlyingType(Multilooking,"logical")} 
        ProcessingStep  {mustBeText,mustBeNonempty,mustBeMember(ProcessingStep, ...
            {'All','Range Chirp','Range Compression','Doppler Centriod','RCMC','Azimuth Compression','Azimuth IFFT and Look Summation'})}
    end
    methods
        function image = process(self)
            %PROCESS Processes the radar datacube to form an image
            %   This method is intended to be overwritten by sub-classes
            %   that define a specific way to process a SAR datacube and
            %   produce an image.

            % TODO SDD: Configure RDA based on 'Operation' property
            switch self.ProcessingStep
                case 'Range Chirp'
                    image = calcRangeChirp();
                case 'All' 
                otherwise
                    warning('Unexpected processing step. No steps completed.')
            end
        end
        function chirpSignal = calcRangeChirp(rda_alg_obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            fs_hz = rda_alg_obj.Radar.Fs_hz;
            bw_hz = rda_alg_obj.Radar.Bandwidth_hz;
            n_samples = rda_alg_obj.Radar.N_pulses;
            chirpSignal = create_lfm_pulse_samples(Causal=0, ...
                                                   ChirpUpDown=1, ...
                                                   Fs_hz=fs_hz, ...
                                                   SweepBandwidth_hz=bw_hz, ...
                                                   NumberOfSamples=n_samples);
        end
        function outputArg = performRangeCompresssion(rda_alg_obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = rda_alg_obj;
        end
        function outputArg = performAzimuthFFt(rda_alg_obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = rda_alg_obj;
        end
        function outputArg = performDopplerCentroidEst(rda_alg_obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = rda_alg_obj;
        end
        function outputArg = performRcmc(rda_alg_obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = rda_alg_obj;
        end
        function outputArg = performAzimuthCompression(rda_alg_obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = rda_alg_obj;
        end
        function outputArg = performAzimuthIfftLookSummation(rda_alg_obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = rda_alg_obj;
        end
        
    end
end

