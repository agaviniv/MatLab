% Operador de Relaxaco

clear all

[x y z]=mangqbts(2);

% Base de Transicao
t{1}=eye(4);
t{2}=2*z{1};
t{3}=2*z{2};
t{4}=4*z{1}*z{2};
t{5}=4*(x{1}*x{2}+y{1}*y{2});
t{6}=4*(x{1}*y{2}-y{1}*x{2});
t{7}=2*x{1};
t{8}=2*y{1};
t{9}=2*x{2};
t{10}=2*y{2};
t{11}=4*x{1}*z{2};
t{12}=4*y{1}*z{2};
t{13}=4*z{1}*x{2};
t{14}=4*z{1}*y{2};
t{15}=4*(x{1}*x{2}-y{1}*y{2});
t{16}=4*(x{1}*y{2}+y{1}*x{2});


%%
% Sistema nuclear
hi = pi*(162*z{2}+(6/2)*z{1}*z{2});
roeq = z{1}+z{2};

%%
gama1=0.0788;
p1=5/6;
gama2=0.0786;
p2=2/3;
E{1} = kron(sqrt(p1)*diag([1 sqrt(1-gama1)]),eye(2));
E{2} = kron(sqrt(p1)*[0 sqrt(gama1); 0 0],eye(2));
E{3} = kron(sqrt(1-p1)*diag([sqrt(1-gama1) 1]),eye(2));
E{4} = kron(sqrt(1-p1)*[0 0;sqrt(gama1) 0],eye(2));
E{5} = kron(eye(2),sqrt(p2)*diag([1 sqrt(1-gama2)]));
E{6} = kron(eye(2),sqrt(p2)*[0 sqrt(gama2); 0 0]);
E{7} = kron(eye(2),sqrt(1-p2)*diag([sqrt(1-gama2) 1]));
E{8} = kron(eye(2),sqrt(1-p2)*[0 0;sqrt(gama2) 0]);

qi1 = 0.8*pi;
qi2 = 0.8*pi;
t1 = 0.2; np = 4; vrod = zeros(16,np);
roin = diag([3 -1 -1 -1]); Dro=zeros(4,4);
for k=1:np
    tk = 2^(k-1)*t1;
    alfa1=1-cos(qi1*tk)^2;
    alfa2=1-cos(qi2*tk)^2;
    F{1} = kron(diag([1 sqrt(1-alfa1)]),eye(2));
    F{2} = kron([0 0;0 sqrt(alfa1)],eye(2));
    F{3} = kron(eye(2),diag([1 sqrt(1-alfa2)]));
    F{4} = kron(eye(2),[0 0;0 sqrt(alfa2)]);
    
    u = expm(-i*hi*tk);
    roin = u*roin*u';
    for m=1:2
        Dro = Dro + F{m}*roin*F{m}';
    end
    Dro2=zeros(4,4);
    for m=1:4
        Dro2 = Dro2 + E{m}*Dro*E{m}';
    end
    Dro3=zeros(4,4);
    for m=3:4
        Dro3 = Dro3 + F{m}*Dro2*F{m}';
    end
    Dro4=zeros(4,4);
    for m=5:8
        Dro4 = Dro4 + E{m}*Dro3*E{m}';
    end       
    roin=Dro4; Dro=zeros(4,4); 
    
    rod = roin - roeq;
    trod=zeros(4,4);
    for kk=1:16
        cotb=0.5*trace(t{kk}*rod);
        trod=trod+cotb*t{kk};
    end
    vrod(:,k) = reshape(trod,16,1);
    
end
