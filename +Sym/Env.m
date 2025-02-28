classdef Env
    properties
        turb
        turbswitch
        ctwind
        g
        rho
    end    
    
    methods
        function obj = Env(Turbulence,wind)
            obj.ctwind=wind;
            obj.g=9.81;
            obj.rho=1.225;
            obj.turb=Aux.dryden([50 50 50],[0.5 0.5 0.5]);
            if Turbulence=="on"
               obj.turbswitch=true;
            else
                obj.turbswitch=false;
            end
        end

        function [Vws,Vwg]=wind(obj,V,dt,euler)
            Vws=Rot.DCMb_ned(euler)'*obj.ctwind';
            Vwg=[0;0;0];
            if obj.turbswitch
                Vwg=obj.turb.gen_turb(V,dt);
            end
        end
    end
end

