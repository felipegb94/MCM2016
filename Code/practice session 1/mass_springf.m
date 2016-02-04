function f = mass_springf(t,y)
f = zeros(4,1);
alpha = 1;
beta = 1;
gamma = 1;
omega = 500000;
f(1) = y(2);
f(2) = alpha*(y(3)-2*y(1)); %acceleration of mass 1
f(3) = y(4);
f(4) = beta*(y(1)-y(3)) + gamma*4*sin(omega*t); %acceleration of mass 2
end