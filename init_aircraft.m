%NomeVariavel=Aircraft.Surface('nome',[x,y,z],c,b,S)
%NomeVariavel.from_vsp('arquivo_polar','arquivo_stab')

wing=Aircraft.Surface('mainwing',[0;0;0],0.482,2.5,1.207);
wing.from_vsp('VSPdata/DV24/WBpolar2.txt','VSPdata/DV24/WBstab2.txt')

ht=Aircraft.Surface('hstab',[1.222;0;0.195],0.259,0.968,0.251);
ht.from_vsp('VSPdata/DV24/Htpolar2.txt','VSPdata/DV24/Htstab2.txt')

vt=Aircraft.Surface('vstab',[0.883;0;0.140],0.262,0.240,0.063);
vt.from_vsp('VSPdata/DV24/Vtpolar3.txt','VSPdata/DV24/Vtstab2.txt')

A=Aero.Aeromod(wing,ht,vt,[0.1;0;0]);
Massmod=Aircraft.Mass(15,[0.1;0;0],diag([0.442 0.657 0.447]));
TF=Aircraft.Thrust(0.075,28.16);

DV=Aircraft.Plane(A,TF,Massmod);

clearvars wing ht vt A Massmod TF