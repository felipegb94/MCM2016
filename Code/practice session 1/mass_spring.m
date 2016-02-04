% simulate two masses on springs

% initial conditions
y0 = [0 0 0 0]';

[t,y] = ode45(@mass_springf,[0,10],y0);

figure(1)
subplot(2,1,1)
plot(t,y(:,1))
title('Position of Mass 1')
xlabel('t*')
ylabel('s1')
subplot(2,1,2)
plot(t,y(:,3))
title('Position of Mass 2')
xlabel('t*')
ylabel('s2')