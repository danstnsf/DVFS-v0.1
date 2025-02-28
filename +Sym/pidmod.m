classdef pidmod<handle
    
    properties
        count
    end
    
    methods
        function obj = pidmod()
            obj.count=0;
        end

        function [u,eo]=exec(obj,xi,xr,ui,xe,xdot)
            [a,b,Va]=Aero.V2a_b(xi(1:3));
            fpa=xi(8)-a*pi/180;
            xie=[0;xi(12);0;0];
            ei=xr-xie;
            %ei(2)*180/pi
            %ei(1)=30*pi/180-xi(9);
            %Controle de guiamento
            kppsi=0.4;%0.2;   %Ganho proporcional da direcao
            kdpsi=0.75;%0.6;   %Ganho derivativo da direcao
            kipsi=0;%0.0001;  %Ganho integratico da direcao
            
            kpphi=0.33;%0.003; %Ganho proporcional do angulo de rolagem
            kdphi=0.16;%0.004; %Ganho derivativo do angulo de rolagem
            kiphi=0.05;     %Ganho integrativo do angulo de rolagem

            %phir=ei(1)*kppsi-xdot(9)*kdpsi+xe(2)*kipsi; %Calculo da referencia do angulo de rolagem adicional
            phir=xr(1)-xi(7);
%             phir=(1+0)*pi/180-xi(7);
%             if obj.count<1
%                 obj.count=obj.count+0.01;
%                 phir=0-xi(7);
%             end
            da=phir*kpphi-xi(4)*kdphi+kiphi*xe(4); %Comando de aileron adicional
           
            %Controle do angulo de derrapagem
            
            kpb=0.0001;%0.001;    %Ganho proporcional do angulo de derrapagem
            kib=0;%0.00001; %Ganho integrativo do angulo de derrapagem
            
            ber=-b;
            dr=ber*kpb;%+xe(3)*kib;
            dr=0;
            %Controle de angulo de arfagem
            kpth=3; %Ganho proporcional
            kdth=0.5;%Ganho derivativo
            
            %Controle de altitude por erro
            kph=-0.02; %Ganho proporcional
            kdh=0.0001;
            kih=-0.00001;%Ganho integrativo

            %Controle de altitude por angulo de voo
            kpfpa=0;%0.5;
            kdfpa=0.1;%0.01

            %thr=xr(2)-xi(8);
            %thr*180/pi
            %thr*180/pi
            thr=(1+0.4309)*pi/180-xi(8);
            de=kpth*thr-kdth*xdot(8);
 
            %de=0;

            %Controle de velocidade
            kpv=0.00005;%0.01;   %Ganho proporcional
            kiv=0; %Ganho integrativo
                
            Vr=16;
            ver=Vr-Va;
            dt=ver*kpv+xe(1)*kiv;
            %dt=0;
            
            a=max(-20*pi/180,min(20*pi/180,ui(1)+da));
            e=max(-30*pi/180,min(30*pi/180,ui(2)-de));
            r=max(-30*pi/180,min(30*pi/180,ui(3)+dr));
            t=max([0 ui(4)+dt]);

            dt=(16-Va)*0.02;

    
            
            de=0;
            u=[da;-0.106444815631119-de;r;1.00612327818113+dt];
            
            eo=[ver;ei(1);ei(2);phir;zeros(8,1)];%[ei(9);ei(7);ber;-ei(12);-ver];
            
        end
    end
end

