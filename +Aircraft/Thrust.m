classdef Thrust
    %THRUST Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Shel
        kmotor

    end
    
    methods
        function obj = Thrust(Shel,kmotor)
            obj.Shel=Shel;
            obj.kmotor=kmotor;
        end
        
        function tvec = thrustvec(obj,rho,V,d)
%             if d<0%||d>1
%                 tvec=zeros(3,1);
%                 return
%             end
            [~,~,Vr]=Aero.V2a_b(V);
            F=0.5*rho*obj.Shel*((obj.kmotor*d)^2-Vr^2)*d/abs(d);
%             if F<0
%                 F=0;
%             end
            tvec=[F;0;0];
        end
    end
end

