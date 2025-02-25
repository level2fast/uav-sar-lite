classdef RadarSat1DataProcessor < SarDataProcessorIf
    %RADAR self class is used to process RadarSat1 SAR data.
    %   Performs digital processing of synthetic apeture radar data.
    %   Defines commonly used methods and properties associated with
    %   general SAR processing and can be configured to operate on a given
    %   radar with a specific data processing algorithm. 
    %   Uses the strategy pattern which enables the use of different 
    %   variants of an algorithm within an object to be able to switch from
    %   one algorithm to another during runtime.
    properties(SetAccess=public,GetAccess=public)
        ImgFormationAlgorithm {mustBeUnderlyingType(ImgFormationAlgorithm,"SarAlgorithmIf")}
        Image     (:,:) 
    end
    properties (Access = protected)
        % Data Processor Params
        Radar           {mustBeUnderlyingType(Radar,"Radar")}
        DataCube  (:,:) {mustBeNonempty} = 0
        Mode            {mustBeMember(Mode,{'stripmap','spotlight'})} = 'stripmap'
        Waveform        {mustBeNonempty} = 0
    end
    methods(Access=public)
        function self = RadarSat1DataProcessor(Algorithm,Radar,DataCube,Mode,Waveform)
            %RADAR Creates a SAR data Processor object 
            self.Radar     = Radar;
            self.DataCube  = DataCube;
            self.Mode      = Mode;
            self.Waveform  = Waveform;
            self.ImgFormationAlgorithm = Algorithm;
        end

        function result = run_img_formation_alg(self)
            %PROCESS Process the radar data using the specified Radar and
            % Algorithm

            % determine strategy that will be used for processing data
                % process data using selected strategy
            try 
               self.Image = self.ImgFormationAlgorithm.process();
               result = true;
            catch 
                error("Failed to process radar data")
            end
        end  
    end
    methods
        function set.Radar(self,radar)
            self.Radar = radar;
        end
        function set.ImgFormationAlgorithm(self,algorithm)
            try mustBeUnderlyingType(algorithm,"SarAlgorithm")
                self.ImgFormationAlgorithm = algorithm;
            catch
                error("ImgFormationAlgorithm must be underlying type SarAlgorithmIf")
            end
        end
    end
end