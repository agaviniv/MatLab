function yprime = ebloch(t,y,gama,B0,B1,w,w3)
% y(1) = q
% y(2) = p
% Parametros

%gama = 2*pi* 100.568/9.39; %carbono MHz/Tesla
% gama = 2*pi*399.952/9.39; %hidrogenio

%B1 = 1e-3; %Tesla
%B0 = 9.39; %Tesla
%w = 1.038*gama*B0;  % sintonia da RF - MHz
%w3 = gama*B0 - w;  % gama*B0 = 100 MHz


% Pulsos nas direcoes X e Y
D = 1;    % +1:X  &  -1:-X
PX_q = gama*B1*D * y(2) * (cos(y(1))/sqrt(1 - y(2)^2));
PX_p = -gama*B1*D * sqrt(1-y(2)^2) * sin(y(1));

D = 0;   % +1:Y  &  -1:-Y
PY_q = gama*B1*D * y(2) * (sin(y(1))/sqrt(1 - y(2)^2));
PY_p = gama*B1*D * sqrt(1-y(2)^2) * cos(y(1));


% Integracao
yprime = [PX_q + PY_q - w3; PX_p + PY_p];

