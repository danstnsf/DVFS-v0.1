function M = DCMa_b(alpha,beta)
    r=0;
    p=(180+alpha)*pi/180;
    h=beta*pi/180;
    M=[cos(p)*cos(h) cos(p)*sin(h) -sin(p);
        -cos(r)*sin(h)+sin(r)*sin(p)*cos(h) cos(r)*cos(h)+sin(r)*sin(p)*sin(h) sin(r)*cos(p);
        sin(r)*sin(h)+cos(r)*sin(p)*cos(h) -sin(r)*cos(h)+cos(r)*sin(p)*sin(h) cos(r)*cos(p)];
end

