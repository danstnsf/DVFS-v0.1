classdef Trim
    properties
        trim_cond="RetoNivelado"
        xtrim
        utrim
        const
        fcost
        A
        B
        Ar
        Br
    end
    methods
        
        function obj=Aerotrim(obj,Fxdot,x0,u0,cts)
            kmax=100;
            fmin=1e-5;
            k=1;
            f=9999;
            opt=optimset('MaxFunEvals',12*500);
            while k<kmax && f>fmin
                if k==1
                    Z=[x0(1:9);u0];
                end
                [Ztrim,f]=fminsearch(@(z) obj.costfun(z,1,cts,Fxdot),Z,opt);
        
                Z=Ztrim;
                k=k+1;
            end
            obj.xtrim=[Ztrim(1:9);x0(10:12)];
            obj.utrim=Ztrim(10:13);
            obj.fcost=f;
            obj.const=cts;
        end

        function f0=costfun(obj,z,k,cts,Fxdot)
            switch k
                case 1 %Straight level flight | cts=[V]
                    x=z(1:9);
                    u=z(10:13);

                    xdot=Fxdot(0,x,u);
                    theta=x(8);
                    V=sqrt(x(1)^2+x(2)^2+x(3)^2);
                    alpha=atan2(x(3),x(1));
                    gam=theta-alpha;
                    
                    Q=[xdot(1:9);V-cts(1);gam;x(2);x(7);x(9)];
                    H=diag(ones(1,14));
                    f0=Q'*H*Q;
            end
            
        end

        function obj=implicitlin(obj,F,xdoto,xo,uo,dxdot,dx,du)
            n=length(xdoto);
            m=length(uo);

            E=zeros(n,n);
            Ap=zeros(n,n);
            Bp=zeros(n,m);
             
            for j=1:n
                xdotplus=xdoto;
                xdotmin=xdoto;

                xdotplus(j)=xdotplus(j)+dxdot;
                xdotmin(j)=xdotmin(j)-dxdot;
                E(:,j)=(F(xdotplus,xo,uo)-F(xdotmin,xo,uo))/(2*dxdot);
            end

            for j=1:n
                xplus=xo;
                xmin=xo;

                xplus(j)=xplus(j)+dx;
                xmin(j)=xmin(j)-dx;
                Ap(:,j)=(F(xdoto,xplus,uo)-F(xdoto,xmin,uo))/(2*dx);
            end

            for j=1:m
                uplus=uo;
                umin=uo;
         
                uplus(j)=uplus(j)+du;
                umin(j)=umin(j)-du;
               
                Bp(:,j)=(F(xdoto,xo,uplus)-F(xdoto,xo,umin))/(2*du);
            end
            
            A=-inv(E)*Ap;
            B=-inv(E)*Bp;

            P=zeros(12,12);
            order=[1 3 5 8 12 2 4 6 7 9 10 11];
            for i=1:12
                P(i,order(i))=1;
            end
            Tr=inv(P);
            Ar=Tr'*A*Tr;
            Br=Tr'*B;

            obj.A=A;
            obj.B=B;
            obj.Ar=Ar;
            obj.Br=Br;

        end


    end
end

