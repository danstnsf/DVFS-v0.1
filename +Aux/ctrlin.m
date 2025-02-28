function [u,time] = ctrlin(time,uo,type,initial_time,amplitude,duration)
    %Control entry types: Aileron Doublet(1), Elevator pulse(2)
    if type==1
        u=repmat(uo,1,length(time));
        return
    end
    [~, id] = min(abs(time - initial_time));
    lenu=length(time);
    ts=time(2)-time(1);
    dursteps=(duration-mod(duration,ts))/ts;
    
    u1=repmat(uo,1,id);
    u3=repmat(uo,1,lenu-id-dursteps);
    switch type
        case 2
            n1=(dursteps-mod(dursteps,2))/2;
            n2=(dursteps+mod(dursteps,2))/2;
            u21=repmat(uo,1,n1)+repmat([amplitude;0;0;0],1,n1);
            u22=repmat(uo,1,n2)+repmat([-amplitude;0;0;0],1,n2);          
            u2=[u21 u22];
        case 3
            u2=repmat(uo,1,dursteps)+repmat([0;amplitude;0;0],1,dursteps);
        case 4
            n1=(dursteps-mod(dursteps,2))/2;
            n2=(dursteps+mod(dursteps,2))/2;
            u21=repmat(uo,1,n1)+repmat([0;amplitude;0;0],1,n1).*linspace(0,1,n1);
            u22=repmat(uo,1,n2)+repmat([0;amplitude;0;0],1,n2).*linspace(1,0,n2);
            u2=[u21 u22];
            
    end
    u=[u1 u2 u3];
end

