function [W espectro] = signw(roin,H0,sw,np,tof,rp,lp)

global Ix Iy Iz sfrq;

dt = 1/(sw);

for tt = 1:np
    U = expm(-i*2*pi*tof*Iz*tt*dt)*expm(-i*H0*tt*dt);
    %U = expm(-i*H0*tt*dt);
    fid(tt) = trace(Iy*(U*roin*U'));
end

espectro = fftshift(fft(fid));
espectro = espectro(1:end/2);
%W = (sw)*(-np/2:(np/2)-1)/np;
W = (sw)*(-np/2:2:(np/2)-1)/np;
