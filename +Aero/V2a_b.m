function [a,b,V]= V2a_b(Vel)
    V=sqrt(Vel(1)^2+Vel(2)^2+Vel(3)^2);
    a=atan2(Vel(3),Vel(1))*180/pi;
    if V==0
        b=0;
    else
        b=asin(Vel(2)/V)*180/pi;
    end
end