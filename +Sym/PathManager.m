classdef PathManager < handle
    properties
        tau=5
        taulead=2;
        Mdp=1.5
        gmax=30*pi/180
        phimax=60*pi/180
        amax=9.81*tan(60*pi/180)
        fpamax=10*pi/180;
        wpswerror=30;

        wpdata
        wpi
        wplen
    end

    methods
        function obj=PathManager(wpfile)
            wpdata=readtable(wpfile);
            obj.wpdata=table2array(wpdata);
            obj.wpi=1;
            obj.wplen=size(obj.wpdata,1);
        end
        
        function xr=get_ref(obj,x,xdot)            
            %Va=x(1:3);va=sqrt(Va'*Va);
            Va=xdot(10:12);va=sqrt(Va'*Va);
            Vg=[1 0 0;0 1 0;0 0 0]*Va;vg=sqrt(Vg'*Vg);
            
            pf=1; %fator de perturbacao
            R=(va*pf)^2/obj.amax; %Raio de curvatura do movimento baseado na velocidade em relacao ao ar + fator de perturbacao
            
            i=obj.wpi; %pega o valor do waypoint atual da memoria do objeto
            pb=x(10:12); %posicao em relacao ao referencial inercial
            wp=[obj.wpdata(:,2:4);obj.wpdata([1 2],2:4)]; %matriz de coordenadas dos waypoints, com repeticao das duas primeiras linhas ao final
            po=wp(i,:)'; %ponto atual (i)
            p1=wp(i+1,:)'; %ponto i+1
            p2=wp(i+2,:)'; %ponto i+2
            
            %Vetores direcao da perna atual e da proxima
            T1=(p1-po);T1=T1/sqrt(T1'*T1);
            T2=(p2-p1);T2=T2/sqrt(T2'*T2);
            
            Tg1=[1 0 0;0 1 0;0 0 0]*T1;
            Tg2=[1 0 0;0 1 0;0 0 0]*T2;
            
            ga=(pi-acos(Tg1'*Tg2))/2;
            p=R/tan(ga)+obj.taulead*vg;
            psw=[1 0 0;0 1 0;0 0 0]*p1-p*Tg1;
            
            perror=sqrt((psw-[1 0 0;0 1 0;0 0 0]*pb)'*(psw-[1 0 0;0 1 0;0 0 0]*pb));
            if perror<obj.wpswerror && i<obj.wplen
                obj.wpi=obj.wpi+1;
                %keyboard
            end
            
            [phir,thetar,en,eh]=obj.getroll_stltrack(po,p1,x,xdot);
            xr=[phir;thetar;en;eh];
        end
        function [phir,thetar,en,eh]=getroll_stltrack(obj,Po,Pi,x,xdot)
            Pb=x(10:12); %posicao atual (NED)
            T=(Pi-Po);T=T/sqrt(T'*T); %vetor direcao da perna atual
            N=[0 -1 0;1 0 0;0 0 0]*T;N=N/sqrt(N'*N);
            B=cross(T,N);
        
            en=N'*(Pb-Po); %erro lateral
            %Vi=Rot.DCMb_ned(x(7:9))*x(1:3);Vg=[1 0 0;0 1 0;0 0 0]*Vi;
            Vi=xdot(10:12);Vg=[1 0 0;0 1 0;0 0 0]*Vi;
            L2=sqrt(Vg'*Vg)*obj.tau;
           
            Ddt=en/tan(obj.gmax);
        
            Ddtmin=min(Ddt,obj.Mdp*L2);
        
            if en<=L2
                Ddtmin=max(Ddtmin,sqrt(L2^2-en^2));
            end
        
            Dwpi=T'*(Pi-Pb);
            Da=Dwpi-Ddtmin;
            Pa=Pi-T*max(0,Da);
            
            Ldir=[1 0 0;0 1 0;0 0 0]*(Pa-Pb);
            c1=cross(Vg,Ldir);
            %c1=cross(Vg,Pa-Pb);
            sineta=sqrt(c1'*c1)/(sqrt(Vg'*Vg)*sqrt(Ldir'*Ldir))*sign(c1(3));
            
            %sineta=sqrt(c1'*c1)/(sqrt(Vg'*Vg)*sqrt((Pa-Pb)'*(Pa-Pb)))*sign(c1(3));
            
            %acmd=2*sqrt(Vg'*Vg)/obj.tau*sineta;
            %roll=max(-obj.phimax,min(atan(acmd/9.81),obj.phimax));
            %keyboard
            ndir=[1 0 0;0 1 0;0 0 0]*T;
            c2=cross(Vg,Ldir);
            sineta2=sqrt(c2'*c2)/(sqrt(Vg'*Vg)*sqrt(Ldir'*Ldir))*sign(c2(3));
            %dpsi=asin(sineta2);
            acmd=2*sqrt(Vg'*Vg)/obj.tau*sineta2;
            phir=max(-obj.phimax,min(atan(acmd/9.81),obj.phimax));

            %Llong=(Pa-[1 0 0;0 1 0;0 0 -1]*Pb)-N'*(Pa-[1 0 0;0 1 0;0 0 -1]*Pb);
            Llong=Pa-[1 0 0;0 1 0;0 0 -1]*Pb;
            eh=Pa(3)+Pb(3);
            Ld=T'*Llong;
            fpa=atan2(eh,sqrt(Ld'*Ld));
            %fpa*180/pi
            %fpa=max(-obj.fpamax,min(fpa,obj.fpamax));
            %vp=Rot.DCMb_ned(xdot(7:9))'*xdot(10:12);
            vp=xdot(10:12);
            thetar=fpa+atan2(x(3),x(1));%-atan2(-vp(3),vp(1))+x(7);
            
            %thetar*180/pi
            %hr=max(-2,min(-Pa(3),2));
            %keyboard
            %amax=9.81*tan(60*pi/180);
            %taulead=4;
            %ucmd=abs(Vi);
            %etamax=min(pi/2,taulead*amax/ucmd/2);
        end
    end
end