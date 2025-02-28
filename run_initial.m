tl=length(S.CFG.time);
us=Aux.ctrlin(S.CFG.time,ut,2,1,2*pi/180,4);
[data,datadot,uctrl]=S.runsym(us,"initial");
plotdata(data,uctrl)