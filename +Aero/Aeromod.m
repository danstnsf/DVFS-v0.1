classdef Aeromod
    
    properties
        Surfaces
        cgloc
        cma
        S
        b

    end
    
    methods
        function obj = Aeromod(Mainwing,HS,VS,cgloc)
            obj.Surfaces={Mainwing,HS,VS};
            obj.cgloc=cgloc;
            obj.cma=Mainwing.cbar;
            obj.S=Mainwing.S;
            obj.b=Mainwing.b;
        end

        function [F,M]=state(obj,V,rho,w,ctrl)
            [alpha,beta,Vr]=Aero.V2a_b(V);
            
            dail=ctrl(1);
            delev=ctrl(2);
            drud=ctrl(3);
            
            Wstate=obj.Surfaces{1}.getstate(alpha,beta,Vr,w,dail);
            
            dw=Wstate(3)*obj.Surfaces{1}.S/obj.Surfaces{1}.b^2/pi*(180/pi);
            
            Hstate=(obj.Surfaces{2}.S/obj.S)*diag([1 1 1 obj.Surfaces{2}.b/obj.b obj.Surfaces{2}.cbar/obj.cma obj.Surfaces{2}.b/obj.b])*obj.Surfaces{2}.getstate(alpha-dw,beta,Vr,w,delev);

            Vstate=(obj.Surfaces{3}.S/obj.S)*diag([1 1 1 obj.Surfaces{3}.b/obj.b obj.Surfaces{3}.cbar/obj.cma obj.Surfaces{3}.b/obj.b])*obj.Surfaces{3}.getstate(alpha,beta,Vr,w,drud);
            
            qinf=Aero.qinf(V,rho);
            Rab=Rot.DCMa_b(alpha,beta);

            Fw=Rab*qinf*obj.S*Wstate(1:3,1);
            Fh=Rab*qinf*obj.S*Hstate(1:3,1);
            Fv=Rab*qinf*obj.S*Vstate(1:3,1);
            
            vposw=diag([-1 1 -1])*obj.Surfaces{1}.ca;
            vposh=diag([-1 1 -1])*obj.Surfaces{2}.ca;
            vposv=diag([-1 1 -1])*obj.Surfaces{3}.ca;
            cg=diag([-1 1 -1])*obj.cgloc;
            
            mcor=obj.S*diag([obj.b obj.cma obj.b]);
            Mw=mcor*qinf*Wstate(4:6,1)+cross(vposw-cg,Fw);
            Mh=mcor*qinf*Hstate(4:6,1)+cross(vposh-cg,Fh);
            Mv=mcor*qinf*Vstate(4:6,1)+cross(vposv-cg,Fv);
            
            F=Fw+Fh+Fv;
            M=Mw+Mh+Mv;
            
        end
    end
end

