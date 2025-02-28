function pinf=qinf(Vel,rho)
    V=sqrt(Vel(1)^2+Vel(2)^2+Vel(3)^2);
    pinf=0.5*V^2*rho;
end