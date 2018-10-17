
clear all;

global tn sfrq H0 Ix Iy Iz ro tof tpwr;

tn = 'H1'; 

[x y z] = mangqbts(3);

Ix = (x{1}+x{2}+x{3});
Iy = (y{1}+y{2}+y{3});
Iz = (z{1}+z{2}+z{3});

sfrq = 500e6;   %Hz

cs = [0 0 0];  %Hz
J = [100 0 0];   %J12 J13 J23 - Hz

Hz = -2*pi*((sfrq-cs(1))*z{1} + (sfrq-cs(2))*z{2} + (sfrq-cs(3))*z{3});
Hi = 2*pi*(J(1)*z{1}*z{2} + J(2)*z{1}*z{3} + J(3)*z{2}*z{3});
H0 = Hz + Hi;

ro = Hz;

tof = 0;
tpwr = 63;  %dB -16:63
pw = 10e-6;   %s
at = 5;  %s
np = 512;  %numero de pontos

clear x y z Hz  cs J;