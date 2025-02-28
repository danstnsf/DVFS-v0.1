classdef adrc < handle
    %ADRC Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        T
        r
        h

        bk
        b1
        b2
        count
    end
    
    methods
        function obj = adrc(ts,r,h)
            obj.T=ts;
            obj.r=r;
            obj.h=h;
            obj.bk=[30;300;1000];
            obj.b1=80;
            obj.b2=-200;
            obj.count=0;
        end
        
        function [vd]=trackdiff(obj,vi,vr)
            a=100;
            Atd=[0 1 0;0 0 1;-a^3 -3*a^2 -3*a];
            Btd=[0;0;a^3];
            vd=Atd*vi+Btd*vr;
        end
        function [zd]=stateobs(obj,y,ui,zi,b)
            e=y-zi(1);
            zd=[0 1 0;0 0 1;0 0 0]*zi+obj.bk*e+[0;b;0]*ui;
        end

        function [da,sd]=latctrl(obj,phi,phiref,dai,xe)
            phiref=1*pi/180;
            [vd]=obj.trackdiff(xe(4:6),phiref);
       
            [zd]=obj.stateobs(phi,dai,xe(1:3),obj.b1);
            
            v1=xe(4);v2=xe(5);

            z1=xe(1);z2=xe(2);z3=xe(3);
            
            e1=v1-z1;e2=v2-z2;
            
            k1=10;
            k2=-11;

            da=(k1*e1+k2*e2-z3)/obj.b1;
            %e1*180/pi
            %da*180/pi

            sd=[zd;vd];
            %da=0;
        end
        function [da,sd]=longctrl(obj,theta,thetaref,dai,xe)
            thetaref=1*pi/180;
            [vd]=obj.trackdiff(xe(4:6),thetaref);
            
            [zd]=obj.stateobs(theta,dai,xe(1:3),obj.b2);
            
            v1=xe(4);v2=xe(5);

            z1=xe(1);z2=xe(2);z3=xe(3);
            
            e1=v1-z1;e2=v2-z2;
            
            k1=10;
            k2=0.5;

            da=(k1*e1+k2*e2-z3)/obj.b2;
           

            sd=[zd;vd];
            
        end

        function [u,eo]=exec(obj,xi,xr,ui,xe,xdoti)
           
            [da,eo1]=obj.latctrl(xi(7),xr(1),ui(1),xe(1:6));
            [de,eo2]=obj.longctrl(xi(8),xr(2),ui(2),xe(7:12));
            
            eo=[eo1;eo2];

            Va=sqrt(xi(1:3)'*xi(1:3));
            dt=(16-Va)*0.02;
            %dt=0;
            %dr=0.1*asin(xi(2)/Va);
            %de=-0.106444815631119;

%             if obj.count<0.5
%                 de=0;
%                 obj.count=obj.count+0.01;
%             end
%             da=0;

            u=[da;de-0.106444815631119;ui(3);1.00612327818113+dt];

            
    
        end

    end
end

