% Resolve a Eq de Bloch, def no arq ebloch.m


gama = 2*pi* 100.568/9.39; %carbono MHz/Tesla
% gama = 2*pi*399.952/9.39; %hidrogenio

B1 = 1e-3; %Tesla
B0 = 9.39; %Tesla
w = 0.99996*gama*B0;  % sintonia da RF - MHz
w3 = gama*B0 - w;  % gama*B0 = 100 MHz


pw = 23.4;      % micro_seg -- carbono:23->pi/2
tspan = [0 4*pw];     



y0 = [0.6999*(3*pi/2);0.6999];     % Condicoes iniciais [q0,p0]

[t,y1] = ode45(@ebloch,tspan,y0,'',gama,B0,B1,w,w3);

for l = 1:size(y1(:,1))
    if ((y1(l,1) < 0) | (y1(l,1) >= 2*pi))
        y1(l,1) = mod(real(y1(l,1)), 2*pi);
    end
end

figure(1)

plot(real(y1(:,1)),real(y1(:,2)),'.b')
axis([0 2*pi -1 1])

%Mx = cos(y(:,1));
%My = sin(y(:,2));


%figure(1)
%plot(t,real(y(:,1)),'r',t,real(y(:,2)),'b')



gama = 2*pi* 100.568/9.39; %carbono MHz/Tesla
% gama = 2*pi*399.952/9.39; %hidrogenio

B1 = 1e-3; %Tesla
B0 = 9.39; %Tesla
w = 0.999954*gama*B0;  % sintonia da RF - MHz
w3 = gama*B0 - w;  % gama*B0 = 100 MHz

pw = 23.4;      % micro_seg -- carbono:23->pi/2
tspan = [0 4*pw];     
y0 = [0.6999*(3*pi/2);0.6999];     % Condicoes iniciais [q0,p0]

w = 1*gama*B0;

[t,y2] = ode45(@ebloch,tspan,y0,'',gama,B0,B1,w,w3);

for l = 1:size(y2(:,1))
    if ((y2(l,1) < 0) | (y2(l,1) >= 2*pi))
        y2(l,1) = mod(real(y2(l,1)), 2*pi);
    end
end

figure(2)

plot(real(y2(:,1)),real(y2(:,2)),'.b')
axis([0 2*pi -1 1])
