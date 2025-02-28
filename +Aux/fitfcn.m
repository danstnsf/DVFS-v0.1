function fn = fitfcn(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10)
    fn=@(x,y) p1+p2*x+p3*y+p4*x^2+p5*x*y+p6*y^2+p7*x^3+p8*x^2*y+p9*x*y^2+p10*y^3;
end

