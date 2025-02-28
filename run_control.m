ct=Sym.pidmod();
wpm=Sym.PathManager('WPdata/oito.txt');
[data,datadot,uctrl]=S.runsym(us,"control",ct,wpm);
plotdata(data,uctrl,wpm)