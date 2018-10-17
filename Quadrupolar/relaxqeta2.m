
clear all

%%
% gama1=linspace(0,1,60);
% p1=2/3;
% gama2=linspace(0,1,60);
% p2=5/6;
% 
% s1 = size(gama1,2); s2 = size(gama2,2);
% roin = 0.25*eye(4); Dro=zeros(4,4);
% etag=zeros(s1,s2);
% for k=1:s1
%     for kk=1:s2
%         Dro = atenuamp(p1,gama1(k),'a',roin);       
%         Dro2 = atenuamp(p2,gama2(kk),'b',Dro);
% 
%         etag(k,kk)=trace((Dro2-0.25*eye(4))^2)/trace((0.25*eye(4))^2);        
%         roin=Dro2;
% 
%     end
% end    
% 
% [G1 G2]=meshgrid(gama1,gama2);
% mesh(G1,G2,etag)
% zlabel('\eta'); ylabel('\gamma_2'); xlabel('\gamma_1')

roin = 0.25*eye(4);
for jj=1:2
%%
p1=5/6;
p2=2/3;
T1a=600e-3; 
T1b=6e-3;

np=70; dt=1e-3;
%roin = 0.25*eye(4);
Dro=zeros(4,4);
eta=zeros(1,np);dd=zeros(np,4);
eta(1)=trace((roin-0.25*eye(4))^2)/trace((0.25*eye(4))^2);
for k=1:np
    gama1=1-exp(-k*dt/T1a);
    gama2=1-exp(-k*dt/T1b);
    Dro = atenuamp(p1,gama1,'a',roin);
    Dro2 = atenuamp(p2,gama2,'b',Dro);
        
    eta(k+1)=trace((Dro2-0.25*eye(4))^2)/trace((0.25*eye(4))^2);    
    roin=Dro2; Dro=zeros(4,4); 
%     if (mod(np,10)==0); roin=rotq(pi,0,'ab',roin); end;
end

plot(0:np,eta,'b.-')
xlabel('time (ms)'); ylabel('\eta');


%%
p1=[1 1 0 0];
p2=[1 0 1 0];

Dro=zeros(4,4); Droc=cell(1,4); d=zeros(np+1,4);
etav=zeros(1,np+1);dd=zeros(np,4);
etav(end)=trace((roin-0.25*eye(4))^2)/trace((0.25*eye(4))^2);
d(end,:)=diag(roin);
for k=np:-1:1
    gama1=1-exp(-k*dt/T1a);
    gama2=1-exp(-k*dt/T1b);
    for ll=1:size(p1,2)
        Dro = atenuamp(p1(ll),gama1,'a',roin);
        Droc{ll} = atenuamp(p2(ll),gama2,'b',Dro);
    end   
    Dro2=0.25*(Droc{1}+Droc{2}+Droc{3}+Droc{4}); d(k,:)=diag(Dro);
    etav(k)=etav(k)+trace((Dro2-0.25*eye(4))^2)/trace((0.25*eye(4))^2);    
    roin=Dro2; Dro=zeros(4,4); 
%     if (mod(np,10)==0); roin=rotq(pi,0,'ab',roin); end;
end
hold on
plot(0:np,etav,'r.-')
axis([-1 np+1 -0.05 max(etav)+0.05*max(etav)])
end