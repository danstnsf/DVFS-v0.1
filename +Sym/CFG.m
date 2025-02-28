classdef CFG

    properties
        ts;
        time;
    end
    
    methods
        function obj = CFG(endtime,timestep)
            if nargin==0
                obj.ts=0.01;
                obj.time=0:0.01:10;
            else
                obj.ts=timestep;
                obj.time=0:timestep:endtime;
            end
        end
    end
end

