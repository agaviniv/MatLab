clear all

% ---- sistema ----
w = 2*pi* [125.0008 125.0011]*1e6;  %homonuclear
j = 2*pi* [200];

I{1,1} = 0.5*kron([0 1;1 0],eye(2));
I{2,1} = 0.5*kron(eye(2),[0 1;1 0]);
I{1,2} = 0.5*kron([0 -i;i 0],eye(2));
I{2,2} = 0.5*kron(eye(2),[0 -i;i 0]);
I{1,3} = 0.5*kron([1 0;0 -1],eye(2));
I{2,3} = 0.5*kron(eye(2),[1 0;0 -1]);

h = -w(1)*I{1,3} -w(2)*I{2,3} + j*I{1,3}*I{2,3};  %Zeeman + acoplamento escalar


% ---- controle ----
 w1 = 2*pi*1e4;
 fase = [pi/2];
 wrf = 2*pi*125e6;
 tp = 25e-6; 

dt= 1e-6; nb=round(tp/dt);
U = cell(1,nb);d=zeros(nb,4);
for kk=0:nb
    hc = -w1*(I{1,1}*cos(wrf*kk*dt+fase)-I{1,2}*sin(wrf*kk*dt+fase))...
         -w1*(I{2,1}*cos(wrf*kk*dt+fase)-I{2,2}*sin(wrf*kk*dt+fase));
    U{kk+1} = expm(-i*(h+hc+wrf*(I{1,3}+I{2,3}))*kk*dt);
    u(kk+1,:)=eig(U{kk+1});
end
ro = U{end}*(I{1,3}+I{2,3})*U{end}';


% % ---- aquisicao ----
% ro = I{1,1}+I{2,1};
%     quad = sqrt(1/1)*[1 -i];
%     Im = quad(1)*(I{1,1}+I{2,1})+quad(2)*(I{1,2}+I{2,2});
%     at = 1; 
%     np = 5000; 
%     tof = 000;
%     wobs = 2*pi*125e6 + 2*pi*tof;
% dt=at/np; eixot=dt*[0:np];
% fid = zeros(1,np);
% for kk=0:np
%     U0 = expm(-i*wobs*(I{1,3}+I{2,3})*dt*kk)*expm(-i*h*dt*kk);
%     fid(kk+1) = trace(U0*ro*U0'*Im)*exp(-dt*kk*150);
% end
% subplot(1,2,1);
% plot (eixot,real(fid),eixot,imag(fid),'r'); xlabel('tempo (s)')
% 
% 
% % ---- processamento ----
% nps = 2^10;
% sw = 1/dt; 
% spect = fft(fid,nps);
% eixof = sw*(0:0.5*nps)/nps + tof;
% subplot(1,2,2);
% plot(eixof,real(spect(1:0.5*nps+1)),eixof,imag(spect(1:0.5*nps+1))); 
% xlabel('frequência (Hz)')

