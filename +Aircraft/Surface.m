classdef Surface
    %SURF Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        name
        cbar %mean aerodynamic chord
        b    %wing span
        S    %wing area
        
        Aero%=Aero.Aerodata() % aerodynamic data class
        pos  %wing leading edge position wrt geometric origin
        ca
        orientation
    end
    
    methods
        function obj = Surface(name,pos,cbar,b,S,orientation)
            obj.Aero=Aero.Aerodata();
            obj.name=name;
            obj.cbar=cbar;
            obj.b=b;
            obj.S=S;
            obj.pos=pos;
            if nargin>5
                if orientation=='v'
                obj.ca=[pos(1)+0.25*cbar;pos(2);pos(3)+0.5*b];
                else
                    obj.ca=[pos(1)+0.25*cbar;pos(2);pos(3)];
                end
            else
                obj.orientation="h";
                obj.ca=[pos(1)+0.25*cbar;pos(2);pos(3)];
            end
        end

        function from_vsp(obj,polar,stab)
            obj.Aero.from_vsp_polar(polar)
            obj.Aero.from_vsp_stab(stab)
        end
        
        function st=getstate(obj,alpha,beta,V,w,d)
            st=obj.Aero.aerostate(alpha,beta,w,V,obj.cbar,obj.b,d);
        end
    end
end

