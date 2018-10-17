% Simulacao para sistema de spin 3/2(quadrupolar)
% Cristal liquido de sodio - 3 linhas


clear all;


Ix = 0.5*[0 sqrt(3) 0 0;sqrt(3) 0 2 0;0 2 0 sqrt(3);0 0 sqrt(3) 0];

Iy = 0.5*i*[0 -sqrt(3) 0 0;sqrt(3) 0 -2 0;0 2 0 -sqrt(3);0 0 sqrt(3) 0];

Iz = 0.5*diag([3 1 -1 -3]);


% Pulso Hard

P90xhard = expm(-i*(pi/20)*Ix);
P90yhard = expm(-i*(pi/20)*Iy);
P90_xhard = expm(i*(pi/20)*Ix);
P90_yhard = expm(i*(pi/20)*Iy);


%------------------------------------------------------------------------------
% Direcao x
% Parametros
	wq = 8;
	tp = pi;
	theta = pi/2;
	fas1 = cos(sqrt(1)*theta/2)*exp(2*i*wq*tp);
   
% Elementos da matriz pulso
fas2 = i*sin(sqrt(1)*theta/2)*exp(i*(2*wq*tp));
fas3 = i*sin(sqrt(1)*theta/2)*exp(i*(2*wq*tp));
fas4 = cos(theta/2)*exp(i*wq*tp);
fas5 = i*sin(theta/2)*exp(i*(wq*tp));
fas6 = i*sin(theta/2)*exp(i*(wq*tp));

% Pulsos de RF - Seletivos

P90x01 = [fas1 fas2 0 0; fas3 fas1 0 0; 0 0 1 0;0 0 0 exp(-i*4*wq*tp)];
P90x12 = [exp(-i*wq*tp) 0 0 0; 0 fas4 fas5 0; 0 fas6 fas4 0; 0 0 0 exp(-i*wq*tp)];
P90x23 = [exp(-i*4*wq*tp) 0 0 0; 0 1 0 0; 0 0 fas1 fas2; 0 0 fas3 fas1];

% Direcao -x
% Parametros
	wq = 8;
	tp = pi;
	theta = pi/2;
	fas1 = cos(sqrt(1)*theta/2)*exp(2*i*wq*tp);
   
% Elementos da matriz pulso
fas2 = i*sin(-sqrt(1)*theta/2)*exp(i*(2*wq*tp));
fas3 = i*sin(-sqrt(1)*theta/2)*exp(i*(2*wq*tp));
fas4 = cos(theta/2)*exp(i*wq*tp);
fas5 = i*sin(-theta/2)*exp(i*(wq*tp));
fas6 = i*sin(-theta/2)*exp(i*(wq*tp));

% Pulsos de RF - Seletivos

P90_x01 = [fas1 fas2 0 0; fas3 fas1 0 0; 0 0 1 0;0 0 0 exp(-i*4*wq*tp)];
P90_x12 = [exp(-i*wq*tp) 0 0 0; 0 fas4 fas5 0; 0 fas6 fas4 0; 0 0 0 exp(-i*wq*tp)];
P90_x23 = [exp(-i*4*wq*tp) 0 0 0; 0 1 0 0; 0 0 fas1 fas2; 0 0 fas3 fas1];


%--------------------------------------------------------------------------------
% Direcao y
% Parametros
	wq = 8;
	tp = pi;
	theta = pi/2;
	fas1 = cos(sqrt(1)*theta/2)*exp(2*i*wq*tp);
   
% Elementos da matriz pulso
fas2 = i*sin(sqrt(1)*theta/2)*exp(i*(2*wq*tp - (pi/2)));
fas3 = i*sin(sqrt(1)*theta/2)*exp(i*(2*wq*tp + (pi/2)));
fas4 = cos(theta/2)*exp(i*wq*tp);
fas5 = i*sin(theta/2)*exp(i*(wq*tp - (pi/2)));
fas6 = i*sin(theta/2)*exp(i*(wq*tp + (pi/2)));

% Pulsos de RF - Seletivos

P90y01 = [fas1 fas2 0 0; fas3 fas1 0 0; 0 0 1 0;0 0 0 exp(-i*4*wq*tp)];
P90y12 = [exp(-i*wq*tp) 0 0 0; 0 fas4 fas5 0; 0 fas6 fas4 0; 0 0 0 exp(-i*wq*tp)];
P90y23 = [exp(-i*4*wq*tp) 0 0 0; 0 1 0 0; 0 0 fas1 fas2; 0 0 fas3 fas1];


% Direcao -y
% Parametros
	wq = 8;
	tp = pi;
	theta = pi/2;
	fas1 = cos(sqrt(1)*theta/2)*exp(2*i*wq*tp);
   
% Elementos da matriz pulso
fas2 = i*sin(-sqrt(1)*theta/2)*exp(i*(2*wq*tp - (pi/2)));
fas3 = i*sin(-sqrt(1)*theta/2)*exp(i*(2*wq*tp + (pi/2)));
fas4 = cos(theta/2)*exp(i*wq*tp);
fas5 = i*sin(-theta/2)*exp(i*(wq*tp - (pi/2)));
fas6 = i*sin(-theta/2)*exp(i*(wq*tp + (pi/2)));

% Pulsos de RF - Seletivos

P90_y01 = [fas1 fas2 0 0; fas3 fas1 0 0; 0 0 1 0;0 0 0 exp(-i*4*wq*tp)];
P90_y12 = [exp(-i*wq*tp) 0 0 0; 0 fas4 fas5 0; 0 fas6 fas4 0; 0 0 0 exp(-i*wq*tp)];
P90_y23 = [exp(-i*4*wq*tp) 0 0 0; 0 1 0 0; 0 0 fas1 fas2; 0 0 fas3 fas1];


% Hamiltoniano
wl = 100;   % Mhz
wq = 0.15;
wos = 150;
H0 = (wos - wl)*Iz;   % referencial girante
Hq = -(1/6)*wq*diag([3 -3 -3 3]);
El90 = expm(-i*(pi/2)*Iz);

% Vetor tempo
sw = 0.1;    %janela de frequencia - MHz
dt = 1/sw;   %incremento no tempo 
np = 2048;   %numero de pontos do FID
t = [0:dt:np*dt];

% Matriz Densidade de Equilibrio
ro = H0+Hq;

% Experimento
Seq1 = P90y12^2;
Seq2 = P90_x01*P90x12^2;
roe = Seq1* ro *Seq1';

% Aplicação de Pulso Hard para Leitura
ro_t = P90_yhard * roe * P90_yhard'; 

for k = 1:1:np;
   O = expm(-i*(H0 + Hq)*t(k));
   FIDp(k) = trace(O*ro_t*O'*(Ix-i*Iy) * exp(-150e-5*t(k)));
   FIDn(k) = trace(O*ro_t*O'*(Ix-i*Iy) * exp(-150e-5*t(k)));
end  

espectro = real(fft(FIDp));
freq = 2*pi*sw*[-(np/2):(np/2)-1]/np;


figure(1)
plot(real(FIDp(:)))
figure(2)
plot(freq,real(espectro(:)))