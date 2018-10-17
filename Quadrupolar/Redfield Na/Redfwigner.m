% A partir das expressoes da matriz densidade de Redfield (artigo Ruben), 
%  escrever a Wigner para cada um dos pontos do espaço de fases

clear all

%Redfield funciona para matriz de desvio, i.e., na aproximacao de altas
%temperaturas!!!!
        %ro = 0.5*diag([3 -1 -1 -1]);
ro = diag([0 0 0 1]); ro=[1;-i;-1;i]*[1 i -1 -i]/4;%ro=[0;1;1;0]*[0 1 1 0]/2;
        %ro = .25*kron([1 1;1 -1],[1 1;1 -1])*ro*kron([1 1;1 -1],[1 1;1 -1])';
        %roeq = .5*diag([3 1 -1 -3]);
%ro=eye(4)/4; ro=[0;1;-1;0]*[0 1 -1 0]/2;
roeq = diag([1/2 1/3 1/6 0]);

p = 1/2;                    %pii = diag(roeq)./(2*diag(roeq)); p=max(pii);

st = 0.5e-3;
C = 1.2e10; J0 = 14.3e-9;  J1 = 4e-9;  J2 = 3.5e-9;

np=36; 

dro = temporo(ro,roeq,np,st,C,J0,J1,J2,p);

for m=1:np;
    dd(m,:)=diag(dro{m});
    eta(m)=4*sum((dd(m,:)-ones(1,4)/4).^2);
    Ml(m)=trace(dro{m}*opmmag(3/2,'z'));
end;

% figure;for m=1:np;subplot(6,6,m);figwig(dro{m});end;
% 
% figure;plot(0:np-1,eta,'o'); ylim([-0.05 0.65]); 
% xlabel('time(ms)'); ylabel('\eta','Rotation',0.0,'FontSize',14);


% % Eta analitico
% vdt=(0:np-1)*st;
% A=(1/4)*sum(diag(ro - roeq)'*[-1; 1; 1; -1])*exp(-2*C*(J1+J2)*vdt);
% B=(1/4)*sum(diag(ro - roeq)'*[1; 1; -1; -1])*exp(-2*C*J2*vdt);
% C=(1/4)*sum(diag(ro - roeq)'*[1; -1; 1; -1])*exp(-2*C*J1*vdt);;
% aneta=16*(A.^2+B.^2+C.^2) - (11/18)*A + (16/3)*B + (8/3)*C + (5/9);
% hold on; plot(1e3*vdt,aneta,1e3*vdt,0.5553*Ml/max(Ml))

%figure;
%prefit(np,0.5,0.5,eta);
% Fitting com canal de atenuacao de amplitude
t = 0:np-1; 
plot(t,eta,'ro'); hold on; h = plot(t,eta,'b'); hold off;
title('Input data'); ylim([0 0.7]);

start=0.23;
options = optimset('Simplex','on','TolX',0.01);
estimated_gama = fminsearch(@(x)disalfafit(x,t,eta,h),start,options)


