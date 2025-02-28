classdef Symod
    %SYMOD Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Plane
        Env
        Trim
        CFG
    end
    
    methods
        function obj = Symod(plane,endt,step,wind,Turb)
            obj.Plane = plane;
            obj.CFG=Sym.CFG(endt,step);
            obj.Env=Sym.Env(Turb,wind);
        end
        
        function obj=trimplane(obj,x0,u0,cts)
           obj.Trim=Sym.Trim;
           Fxd=@(xdot,x,u) obj.xdot(0,x,u,true)-xdot;
           obj.Trim=obj.Trim.Aerotrim(Fxd,x0,u0,cts);
        end
        
        function [data,dotdata,uctrl,error,nerror]=runsym(obj,us,type,controller,manager)
            xo=obj.Trim.xtrim;
            xo(10)=0;
            nint=length(obj.CFG.time);
            dt=obj.CFG.ts;
            xsym=zeros(nint,12);
            dotdata=zeros(nint,12);
            xi=xo;
            u=[us us(:,end)];
            eo=zeros(1,1);
            xe=zeros(12,1);
            xsi=[xi;xe];
            error=zeros(1,nint);
            nerror=zeros(2,nint);
            
            if type=="control"
                for i=1:nint
                    ui=u(:,i);
                    xdoti=obj.xdot(0,xi,ui);
                    xr=manager.get_ref(xi,xdoti);
                    [uc,eo]=controller.exec(xi,xr,ui,xe,xdoti);
                    [~,xs]=ode45(@(t,x) obj.xdots(t,x,uc,eo),[0 dt/2 dt],xsi);
                    xsi=xs(3,:)';
                    xi=xsi(1:12);
                    xsym(i,:)=xi';
                    xe=xsi(13:24);
                    error(1:12,i)=xe;
                    u(:,i:i+1)=[uc uc];
                    nerror(:,i)=xr(3:4)';
                end

            elseif type=="initial"
                 for i=1:nint
                    ui=us(:,i);
                    [t,x]=ode45(@(t,x) obj.xdot(t,x,ui),[0 dt/2 dt],xi);
                    dotdata(i,:)=obj.xdot(0,xi,ui);
                    xsym(i,:)=x(3,:);
                    xi=x(3,:)';
                 end
            else

                 for i=1:nint
                    ui=u(:,i);
                    [t,x]=ode45(@(t,x) obj.xdot(t,x,ui),[0 dt/2 dt],xi);
                    dotdata(i,:)=obj.xdot(0,xi,ui);
                    xsym(i,:)=x(3,:);
                    xi=x(3,:)';
                end
            end
            data=[[obj.CFG.time obj.CFG.time(end)+obj.CFG.ts]' [xo';xsym]];
            uctrl=u;
        end
        function obj=linearize(obj,dxdot,dx,du)
           xo=obj.Trim.xtrim;
           uo=obj.Trim.utrim;
            
           xdoto=zeros(length(xo),1);
           F=@(xdot,x,u) obj.xdot(0,x,u,true)-xdot;
           obj.Trim=obj.Trim.implicitlin(F,xdoto,xo,uo,dxdot,dx,du);
        end
        
        function xdot = xdot(obj,t,x,u,istrim)
            rho=obj.Env.rho;
            g=obj.Env.g;
            m=obj.Plane.Mass.mass;
            I=obj.Plane.Mass.Im;
            w=x(4:6);
            euler=x(7:9);
            
            Va=sqrt(x(1)^2+x(2)^2+x(3)^2);
            V=x(1:3);
            if nargin<5
                [Vws,Vwg]=obj.Env.wind(Va,obj.CFG.ts,euler);
                V=V+Vws+Vwg;
            end

            [Fa,Ma]=obj.Plane.Aero.state(V,rho,w,u);
            Ft=obj.Plane.Thrust.thrustvec(rho,V,u(4));
            Fg=Rot.DCMb_ned(euler)'*[0;0;g]*m;
            Fb=Fa+Ft+Fg;
            Mb=Ma;

            Vdot=Fb/m-cross(w,V);

            Wdot=I\(Mb-cross(w,I*w));

            Eulerdot=Rot.eulerdot(euler)*w;
            
            P_edot=Rot.DCMb_ned(euler)*V;

            xdot=[Vdot;Wdot;Eulerdot;P_edot];
        end
        function xdot=xdots(obj,t,x,u,erri)
            xd1=obj.xdot(t,x,u);
            xdot=[xd1;erri];
        end
    end
end

