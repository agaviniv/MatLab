clear all

% |AB>
% Porta rotacao de angulo theta na direcao y controlada pelo A
ryA = @(theta) [1 0 0 0;0 1 0 0;0 0 cos(theta) -sin(theta)...
                               ;0 0 sin(theta) cos(theta)];

xba = [1 0 0 0;0 0 0 1;0 0 1 0;0 1 0 0];


% Porta rotacao de angulo theta na direcao y controlada pelo A
ryB = @(theta) [1 0 0 0;0 cos(theta) 0 -sin(theta);0 0 1 0;...
                               ;0 sin(theta) 0 cos(theta)];

xab = [1 0 0 0;0 1 0 0;0 0 0 1;0 0 1 0];


% Circuito de AD para A e B
AD=@(t1,t2) xab*ryB(t2)*xba*ryA(t1);
u=[1:4];
for kk=1:4
 pp(1,:)=double(u==kk);
for k=1:128;
    s=k*pi/128;
    pp(k+1,:)=diag(AD(s,s)*diag(pp(1,:))*AD(s,s)');
    eta(k)=sum((pp(k,:)-.25*[1 1 1 1]).^2);
end;
% figure;
% plot(pp,'.')
end;


% Circuito de AD para A e B
ADA=@(t) xba*ryA(t);
ADB=@(t) xab*ryB(t);
u=[1:4];
for kk=1:4
 pp(1,:)=double(u==kk);
for k=1:128;
    sa=k*pi/53; sb=k*pi/7;    
    pp(k+1,:)=diag(kron(trlsq(ADA(sa)*diag(pp(1,:))*ADA(sa)')...
                   ,trmsq(ADB(sb)*diag(pp(1,:))*ADB(sb)')));

end;
% figure;
% plot(pp,'.')
end


% Circuito de PD para A e B
PD=@(t) ryA(t);
for p=0:0
    ket=sqrt(1/4)*[1;exp(i*pi*p/2);exp(i*pi*p);exp(-i*pi*p/2)];
    ro=ket*ket';
    for k=1:8;
        figure;bar3(ro);zlim([-1 1])
        sa=k*pi/8; sb=k*pi/8;    
%         keta = PDA(sa)*ket; roa=keta*keta';
%         ketb = PDB(sb)*keta; rob=ketb*ketb';
%         ro = kron(trlsq(roa),trmsq(rob)); trace(ro)
%         coe(k,:) = [ro(1,2) ro(1,3) ro(1,4) ro(2,3) ro(2,4) ro(3,4)];
        
       roa = trlsq(PD(sa)*kron(trlsq(ro),[1 0;0 0])*PD(sa)');
       rob = trmsq(PD(sb)*kron(trmsq(ro),[1 0;0 0])*PD(sb)');
       ro = kron(roa,rob); 
       coe(k,:) = [ro(1,2) ro(1,3) ro(1,4) ro(2,3) ro(2,4) ro(3,4)];
       
    end;
%     figure;
%     plot(coe,'.')
end;


% Descricao via Wigner do PD
p=0;
ket=sqrt(1/4)*[1;exp(i*pi*p/2);exp(i*pi*p);exp(-i*pi*p/2)];
ro=ket*ket';
for k=1:8;
    figwig(kron(PD(sa)*kron(trlsq(ket*ket'),[1 0;0 0])*PD(sa)',...
            PD(sb)*kron(trmsq(ket*ket'),[1 0;0 0])*PD(sb)'));
    pause(0.8);
end
        