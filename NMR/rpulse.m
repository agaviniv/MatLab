function [matpulso ropulsado] = rpulse(pw,ph,pre_delay,pos_delay)

global tn H0 Ix Iy Iz ro tof tpwr sfrq;

if (exist('tpwr') == 0)
    'Declare o parametro tpwr para a potencia do pulso'
    return;
end
if ((tpwr < -16) | (tpwr > 63))
    'Valor de tpwr deve estar entre -16 e 63'
    return;
end


if ((exist('Ix') == 0) || (exist('Iy') == 0) || (exist('Iz') == 0))
    'Declare Ix, Iy e Iz'
    return;
end

if (exist('tn') == 0)
    'Declare o parametro tn para selecionar o nucleo'
    return;
else
    switch tn
        case 'H1'
            gama = 26.75e7;  %rad/T.s
        case 'C13'
            gama = 6.73e7;
    end
end

if (exist('tof') == 0)
    'Declare o parametro tof para off-set de frequencia do pulso'
    return;
end

%u_pw360=4e-7;
%B1m = 2*pi/(gama*pw360);
Bref = 4e6/gama;  %tpwr=63 referencia: angulo 360; tempo 4us
B = Bref/10^((63-tpwr)/20);

Hrf = -gama*B*(cos(ph)*Ix - sin(ph)*Iy); %on-resonance
%Hrf = -pi/(2*pw)*(cos(ph)*Ix - sin(ph)*Iy);

Utof = expm(i*2*pi*(sfrq+tof)*Iz*pw);
matpulso = expm(-i*Utof*(Hrf+H0)*Utof'*pw);
%U = expm(i*2*pi*tof*pw*Iz)*matpulso*expm(-i*2*pi*tof*pw*Iz);

ropulsado = matpulso * ro * matpulso';
