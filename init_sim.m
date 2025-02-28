
Vr=11;
alpha=0;    
beta=0;

V=[Vr*cosd(alpha)
    Vr*sind(beta)
    Vr*sind(alpha)];

w=[0*pi/180;0*pi/180;0*pi/180];

eul=[0
    0
    0];

pos=[0;0;0];

u=[0*pi/180
    0*pi/180
    0*pi/180
    0];

xo=[V;w;eul;pos];
uo=[0;0;0;0];

clearvars u eul pos w V alpha beta Vr

S=Sym.Symod(DV,20,0.01,[0 0 0],"off");
S=S.trimplane(xo,uo,[10]);
S.Trim.fcost
ut=S.Trim.utrim;
xt=S.Trim.xtrim;

clearvars xo uo