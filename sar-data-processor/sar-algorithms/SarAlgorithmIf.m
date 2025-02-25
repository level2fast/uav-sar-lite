classdef (Abstract) SarAlgorithmIf
    %SARALGORITHM This class defines general properties and mehthods used across
    % SAR algorithms.
    %   This class is an interface that defines general methods used for processing
    %   SAR data and properties that can be used to configure a SAR processor. 
    properties(GetAccess=public,SetAccess=protected)
        ImageData (:,:)
        Name {mustBeMember(Name,{'Range Doppler','Chirp Scaling','SPECAN','OMEGA-K'})} = 'Range Doppler'
    end
    properties(Access=public)
        Mode {mustBeMember(Mode,{'stripmap','spotlight'})} = 'stripmap'
        Radar           {mustBeUnderlyingType(Radar,"Radar")} 
        Datacube  (:,:) {mustBeNonempty} = 0    
    end
    methods(Abstract)
        image = process(sar_alg_obj,radar,datacube)
            %PROCESS Processes the radar datacube to form an image
            %   This method is intended to be overwritten by sub-classes
            %   that define a specific way to process a SAR datacube and
            %   produce an image.
    end
end

