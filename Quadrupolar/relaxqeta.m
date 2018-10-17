% matriz densidade a partir dos dados do Ruben

clear all; %close all

%______________________________________
%%%%%%%%%%%%%  Redfield  %%%%%%%%%%%%%%
Iz = 0.5*diag([3 1 -1 -3]);
I = 3/2;
wq = 16e3;
eta = 0;

C = 1.2e10;

roeq = (2*Iz+3*eye(4))/12;%expm(-Iz);(0.6029*Iz);
% 3.61587 - fator que multiplicado com roeq produz eta experimental
ro = 0.25*eye(4); 
roeq=Iz;ro=[0 .5 .5 .5;.5 0 .5 .5;.5 .5 0 .5;.5 .5 .5 0];

J1=3.8e-9; J2=3.4e-9; 

np=80;
eta=zeros(1,np); dd=zeros(np,4);
eta(1)=trace((ro-0.25*eye(4))^2)/trace((0.25*eye(4))^2);
dt = 1e-3;
for kk = 1:np
    R01 = trace(ro - roeq);
    R02 = sum(diag(ro - roeq)'*[-1; 1; 1; -1]);
    R03 = sum(diag(ro - roeq)'*[1; 1; -1; -1]);
    R04 = sum(diag(ro - roeq)'*[1; -1; 1; -1]);
    Rv = [R01 R02*exp(-2*C*(J1+J2)*kk*dt) R03*exp(-2*C*J2*kk*dt) R04*exp(-2*C*J1*kk*dt)];
    drot(1) = roeq(1,1)+(1/4)*(Rv*[1;-1;1;1]);
    drot(2) = roeq(2,2)+(1/4)*(Rv*[1;1;1;-1]);
    drot(3) = roeq(3,3)+(1/4)*(Rv*[1;1;-1;1]);
    drot(4) = roeq(4,4)+(1/4)*(Rv*[1;-1;-1;-1]);
    dd(kk,:) = drot;
    ro = diag(drot);
        
    eta(kk+1)=trace((ro-0.25*eye(4))^2)/trace((0.25*eye(4))^2);
end    

%_________________________________________________
%%%%%%%%%%%% Atenuação de Amplitude %%%%%%%%%%%%%%
gama1=0.0620;
p1=5/6;
gama2=0.1567;
p2=2/3;
E{1} = kron(sqrt(p1)*diag([1 sqrt(1-gama1)]),eye(2));
E{2} = kron(sqrt(p1)*[0 sqrt(gama1); 0 0],eye(2));
E{3} = kron(sqrt(1-p1)*diag([sqrt(1-gama1) 1]),eye(2));
E{4} = kron(sqrt(1-p1)*[0 0;sqrt(gama1) 0],eye(2));
E{5} = kron(eye(2),sqrt(p2)*diag([1 sqrt(1-gama2)]));
E{6} = kron(eye(2),sqrt(p2)*[0 sqrt(gama2); 0 0]);
E{7} = kron(eye(2),sqrt(1-p2)*diag([sqrt(1-gama2) 1]));
E{8} = kron(eye(2),sqrt(1-p2)*[0 0;sqrt(gama2) 0]);

roin = 0.25*eye(4); Dro=zeros(4,4);
eta2=zeros(1,np);dd2=zeros(np,4);
eta2(1)=trace((roin-0.25*eye(4))^2)/trace((0.25*eye(4))^2);
for k=1:np
    for i=1:4
        Dro = Dro + E{i}*roin*E{i}';
    end
    Dro2=zeros(4,4);
    for i=5:8
        Dro2 = Dro2 + E{i}*Dro*E{i}';
    end    
    eta2(k+1)=trace((Dro2-0.25*eye(4))^2)/trace((0.25*eye(4))^2);
    dd2(k,:)=diag(Dro2);
    roin=Dro2; Dro=zeros(4,4); 
end   

%____________________________________________
%%%%%%%%%%%% Atenuação de Fase %%%%%%%%%%%%%%
qi1 = 0.8*pi;
qi2 = 0.8*pi;
dt = 1e-3;
roin = 0.25*eye(4); Dro=zeros(4,4);
eta3=zeros(1,np);dd3=zeros(np,4);
eta3(1)=trace((roin-0.25*eye(4))^2)/trace((0.25*eye(4))^2);
for k=1:np
    alfa1=1-cos(qi1*dt*k)^2;
    alfa2=1-cos(qi2*dt*k)^2;
    F{1} = kron(diag([1 sqrt(1-alfa1)]),eye(2));
    F{2} = kron([0 0;0 sqrt(alfa1)],eye(2));
    F{3} = kron(eye(2),diag([1 sqrt(1-alfa2)]));
    F{4} = kron(eye(2),[0 0;0 sqrt(alfa2)]);
    for i=1:2
        Dro = Dro + F{i}*roin*F{i}';
    end
    Dro2=zeros(4,4);
    for i=3:4
        Dro2 = Dro2 + F{i}*Dro*F{i}';
    end
    
    eta3(k+1)=trace((Dro2-0.25*eye(4))^2)/trace((0.25*eye(4))^2);
    dd3(k,:)=diag(Dro2);
    roin=Dro2; Dro=zeros(4,4); 
end   

%________________________________________________________
%%%%%%%%%%%% Atenuação de Amplitude e Fase %%%%%%%%%%%%%%
qi1 = 0.8*pi;
qi2 = 0.8*pi;
dt = 1e-3;
roin = 0.25*eye(4); Dro=zeros(4,4);
eta4=zeros(1,np);dd4=zeros(np,4);
eta4(1)=trace((roin-0.25*eye(4))^2)/trace((0.25*eye(4))^2);
for k=1:np
    alfa1=1-cos(qi1*dt*k)^2;
    alfa2=1-cos(qi2*dt*k)^2;
    F{1} = kron(diag([1 sqrt(1-alfa1)]),eye(2));
    F{2} = kron([0 0;0 sqrt(alfa1)],eye(2));
    F{3} = kron(eye(2),diag([1 sqrt(1-alfa2)]));
    F{4} = kron(eye(2),[0 0;0 sqrt(alfa2)]);
    for i=1:2
        Dro = Dro + F{i}*roin*F{i}';
    end
    Dro2=zeros(4,4);
    for i=1:4
        Dro2 = Dro2 + E{i}*Dro*E{i}';
    end
    Dro3=zeros(4,4);
    for i=3:4
        Dro3 = Dro3 + F{i}*Dro2*F{i}';
    end
    Dro4=zeros(4,4);
    for i=5:8
        Dro4 = Dro4 + E{i}*Dro3*E{i}';
    end
    
    eta4(k+1)=trace((Dro2-0.25*eye(4))^2)/trace((0.25*eye(4))^2);
    dd4(k,:)=diag(Dro4);
    roin=Dro4; Dro=zeros(4,4); 
end

eixo=[0:np];figure(1);
plot(eixo,[eta; eta2; eta3; eta4])
axis([0 np -0.2 0.7]); ylabel('\eta'); xlabel('tempo (ms)');


%________________________________________________________
%%%%%%%%%%%%%%  Despolarizacao  %%%%%%%%%%%%%%
p1=0.95; p2=0.95;
roin = 0.25*eye(4); Dro=zeros(4,4); Dro2=zeros(4,4);
eta5=zeros(1,np);dd5=zeros(np,4);
eta5(1)=trace((roin-0.25*eye(4))^2)/trace((0.25*eye(4))^2);
for k=1:np
    
    Dro = despolar(p1,'a',roin);       
    Dro2 = despolar(p2,'b',Dro);
    
    eta5(k+1)=trace((Dro2-0.25*eye(4))^2)/trace((0.25*eye(4))^2);
    dd5(k,:)=diag(Dro2);
    roin=Dro2;
end


%___________________
%%%%%%


