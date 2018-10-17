function [tt fff] = signt(roin,H0,at,np,tof,fasrcv)
% Gerador de sinal da magetizacao no dominio temporal
%
%   roin   : matriz densidade inicial
%   H0     : Hamiltoniana(o) livre
%   at     : tempo de aquisicao do sinal
%   np     : numero de pontos - integracao do sinal
%   tof    : frequencia de "off-set" do receptor
%   fasrcv : fase do receptor

d = log(size(roin,1))/log(2);  %dimensao do espaco de spin

%[Ix Iy Iz] = mangqbts(d);      %operadores de spin
global Ix Iy Iz;


% ===========  SISTEMAS  ===========

% I - homonuclear 2 spin-1/2
%Iobs =  (Ix{1}+Ix{2}) + exp(i*fasrcv)*(Iy{1}+Iy{2});      

% IIa - heteronuclear 2 spin-1/2 - medindo A (1)
%Iobs = exp(i*fasrcv)*Ix{1} + exp(i*fasrcv)*Iy{1};        
%Izof = Iz{1};  

% IIb - heteronuclear 2 spin-1/2 - medindo B (2)
%Iobs = exp(i*fasrcv)*Ix{end} + exp(i*fasrcv)*Iy{end};        
%Izof = Iz{end};

% III - quadrupolar spin-3/2
%Iobs = exp(i*fasrcv)*opmmag(3/2,'x');  
%Izof=opmmag(3/2,'z')     


dt = at/np;

for tt = 1:np
    %U = expm(-i*2*pi*tof*Iz*tt*dt)*expm(-i*H0*tt*dt);
    U = expm(-i*H0*tt*dt);
    fff(tt) = trace(Iy*(U*roin*U'));
end

tt = [0:dt:at-dt];