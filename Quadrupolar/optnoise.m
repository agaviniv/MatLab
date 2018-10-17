clear all
% %______________________________________
% %%%%%%%%%%%%%  Redfield  %%%%%%%%%%%%%%
Iz = 0.5*diag([3 1 -1 -3]);
I = 3/2;
wq = 16e3;
eta = 0;
C = 1.2e10;
roeq = (2*Iz+3*eye(4))/12;
ro = 0.25*eye(4); 

R01 = trace(ro - roeq);
R02 = sum(diag(ro - roeq)'*[-1; 1; 1; -1]);
R03 = sum(diag(ro - roeq)'*[1; 1; -1; -1]);
R04 = sum(diag(ro - roeq)'*[1; -1; 1; -1]);

J1=3.8e-9; J2=3.4e-9; 

np=40;
eta=zeros(1,np);
eta(1)=trace((ro-0.25*eye(4))^2)/trace((0.25*eye(4))^2);
dt = 1e-3;
for kk = 1:np
    eta(kk)=trace((ro-0.25*eye(4))^2)/trace((0.25*eye(4))^2);
    Rv = [R01 R02*exp(-2*C*(J1+J2)*kk*dt) R03*exp(-2*C*J2*kk*dt) R04*exp(-2*C*J1*kk*dt)];
    drot(1) = roeq(1,1)+(1/4)*(Rv*[1;-1;1;1]);
    drot(2) = roeq(2,2)+(1/4)*(Rv*[1;1;1;-1]);
    drot(3) = roeq(3,3)+(1/4)*(Rv*[1;1;-1;1]);
    drot(4) = roeq(4,4)+(1/4)*(Rv*[1;-1;-1;-1]);
   
    ro = diag(drot);
            
end    
eta=eta+rand(1,np)/100;
figure(10);plot([0:np-1],eta,'.-');ylim([0 0.7])
etat = eta; %clear eta;

 etaexp=load('etaexp');
 eta=etaexp.eta;    
 eta=eta*0.5132/6.7531;    eta([2 5 6 7])=etat([2 5 6 7]);%
 np = size(eta,2);

% Fitting com canal de atenuacao de amplitude
t = 0:np-1; 
figure;
plot(t,eta,'ro'); hold on; h = plot(t,eta,'b'); hold off;
title('Input data'); ylim([0 0.7])

p(1)=2/3;
p(2)=5/6;
start = [0.08;0.08]; %gama's
options = optimset('Simplex','on','TolX',0.001);
estimated_gama = fminsearch(@(x)etafit(x,t,eta,p,h),start,options)


% Fitting com redfield
t = 0:np-1; 
figure;
plot(t,eta,'ro'); hold on; h = plot(t,eta,'b'); hold off;
title('Input data'); ylim([0 0.7])


start = [3.8e-9;3.4e-9]; %J's
options = optimset('Simplex','on','TolX',0.001);
estimated_J = fminsearch(@(x)etafitJ(x,t,eta,h),start,options)