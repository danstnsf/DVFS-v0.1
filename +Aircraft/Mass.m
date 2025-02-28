classdef Mass

    properties
        CGloc
        mass
        Im
    end
    
    methods
        function obj = Mass(mass,cg,im)
            obj.CGloc=cg;
            obj.mass=mass;
            obj.Im=im;
        end
    end
end

