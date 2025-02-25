classdef (Abstract) SarDataProcessorIf < handle
    %SarDataProcessor This class is the base class for all SAR data processors.
    properties(Abstract,SetAccess=public,GetAccess=public)
        ImgFormationAlgorithm {mustBeUnderlyingType(ImgFormationAlgorithm,"SarAlgorithmIf")}
        Image     (:,:) 
    end
    methods(Abstract)
        result = run_img_formation_alg(data)
    end
    methods(Access=protected)
        function outputArg = calcGroundResolution(object)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = object;
        end
        function outputArg = calcGroundBasis(object)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = object;
        end
        function outputArg = calcIntegrationTime(object)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = object;
        end
        function outputArg = calcSlantRangeToGroundRange(object)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = object;
        end
        function outputArg = calcCrossRangeResolution(object)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = object;
        end
        function range_resolution = calcRangeResolution(object)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            range_resolution = object;
        end
    end
end