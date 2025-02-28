function edot= eulerdot(euler)
p=euler(1);
t=euler(2);
edot=[1 sin(p)*tan(t) cos(p)*tan(t);0 cos(p) -sin(p);0 sin(p)/cos(t) cos(p)/cos(t)];
end

