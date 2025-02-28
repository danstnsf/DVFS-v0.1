classdef dryden < handle
    properties
        len
        stdev
        hist
        count
    end
    methods
        function obj=dryden(len,stdev)
            obj.len=len;
            obj.stdev=stdev;
            obj.hist=zeros(1000,8);
            obj.count=1;
        end
            
        function Vws = gen_turb(obj,V,dt)
            i=obj.count;
            Dx=V*dt;
            sigma_wn=[obj.stdev(1)*sqrt(2*obj.len(1)/Dx);obj.stdev(2)*sqrt(obj.len(2)/Dx);obj.stdev(3)*sqrt(obj.len(3)/Dx)];
            xin=diag(sigma_wn)*randn(3,1);
            obj.hist(i,4:6)=xin';

            R=Dx*[1/(2*obj.len(1));1/(2*obj.len(2));1/(2*obj.len(3))];
        
            A=[(1-R(1))/(1+R(1));(1-R(2))/(1+R(2));(1-R(3))/(1+R(3))];
            B=[R(1)/(1+R(1));R(2)/(1+R(2));R(3)/(1+R(3))];
        
            u=A(1)*obj.hist(i,1)+2*B(1)*obj.hist(i,4);

            yv=A(2)*obj.hist(i,7)+2*B(2)*obj.hist(i,5);
            yw=A(3)*obj.hist(i,8)+2*B(3)*obj.hist(i,6);

            v=A(2)*obj.hist(i,2)+B(2)*(yv+obj.hist(i,7))+sqrt(3)/(1+R(2))*(yv-obj.hist(i,7));
            w=A(3)*obj.hist(i,3)+B(3)*(yw+obj.hist(i,8))+sqrt(3)/(1+R(3))*(yw-obj.hist(i,8));
            
            obj.hist(i+1,:)=[u v w xin' yv yw];
            
            Vws=[u v w]';
            obj.count=obj.count+1;
        end
    end
end

