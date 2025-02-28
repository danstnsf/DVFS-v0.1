function [rollcmd]=getroll_stltrack(Po,Pi,x)
    Pb=x(10:12);
    T=(Pi-Po);T=T/sqrt(T'*T);
    N=[0 -1 0;1 0 0;0 0 0]*T;N=N/sqrt(N'*N);
    %B=cross(T,N);

    en=N'*(Pb-Po);
    Vi=Rot.DCMb_ned(x(7:9))*x(1:3);Vg=[1 0 0;0 1 0;0 0 0]*Vi;
    
    tau=10;
    L2=sqrt(Vg'*Vg)*tau;
    
    gmax=30*pi/180;
    Ddt=en/tan(gmax);

    Mdp=1.5;
    Ddtmin=min(Ddt,Mdp*L2);

    if en<=L2
        Ddtmin=max(Ddtmin,sqrt(L2^2-en^2));
    end

    Dwpi=T'*(Pi-Pb);
    Da=Dwpi-Ddtmin;
    Pa=Pi-T*max(0,Da);
    
    c1=cross(Vg,Pa-Pb);
    sineta=sqrt(c1'*c1)/(sqrt(Vg'*Vg)*sqrt((Pa-Pb)'*(Pa-Pb)))*sign(c1(3));
    
    acmd=2*sqrt(Vg'*Vg)/tau*sineta;
    rollcmd=max(-60*pi/180,min(atan(acmd/9.81),60*pi/180));

    %amax=9.81*tan(60*pi/180);
    %taulead=4;
    %ucmd=abs(Vi);
    %etamax=min(pi/2,taulead*amax/ucmd/2);
end