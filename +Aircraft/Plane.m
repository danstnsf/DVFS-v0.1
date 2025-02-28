classdef Plane
    %PLANE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Aero
        Thrust
        Mass
    end
    
    methods
        function obj = Plane(aero,thrust,mass)
            obj.Aero=aero;
            obj.Thrust=thrust;
            obj.Mass=mass;
        end
    end
end

