%Funcao de fidelidade
function F = Fid(A,U,S,nf)

global wq ix iy iz iz2 L D f0 f1 dt a1 a2 a3 a4

U1 = eye(2*S+1);
penal = 0;
for s = 1:nf
%Calculo da penalidade:
    if ((A(s,1)<L(1,1)) | (A(s,1)<L(2,1)) | (A(s,3)<L(1,3)) | (A(s,3)<L(2,3)))
        penal = 1;
    end
    Ut = U1;
    
    %Calculo da evolucao unitaria:
    Hrf = -A(s,1)*A(s,3)*(ix*cos(A(s,2)) + iy*sin(A(s,2)));
    Hq = A(s,3)*(wq/2)*iz2;
    U0 = expm(-i*(Hrf + Hq));
    U1 = U0*U1;
    
    %Diferenciais para calculo de erro:
    r = 4*(s-1);
    B(r+1) = A(s,1)*A(s,3) + A(s,1)*D(3) + D(1)*A(s,1)*A(s,3);
    B(r+2) = A(s,2) + D(2);
    B(r+3) = wq*A(s,3) + D(4)*A(s,3) + wq*D(3);
    B(r+4) = A(s,3)*D(5);
    
    %Operadores de desvio:
    G(:,:,r+1) = Ut*(expm(i*B(r+1)*(ix*cos(A(s,2)) + iy*sin(A(s,2))) - i*Hq) - U0)*U1';
    G(:,:,r+2) = Ut*(expm(i*A(s,1)*A(s,3)*(ix*cos(B(r+2)) + iy*sin(B(r+2))) - i*Hq) - U0)*U1';
    G(:,:,r+3) = Ut*(expm(-i*Hrf - i*iz2*B(r+3)/2) - U0)*U1';
    G(:,:,r+4) = Ut*(expm(-i*Hrf - i*Hq + i*B(r+4)*iz) - U0)*U1';
end

%Evolucao e erro durante o tempo morto:
U0 = expm(-i*dt*(wq/2)*iz2);
e3 = U1*(expm(-i*iz2*(wq*dt + D(4)*dt + wq*D(3))/2) - U0);
e4 = U1*(expm(i*dt*D(5)*iz - i*iz2*wq*dt/2) - U0);

U1 = U0*U1;
p = trace(U'*U1);
%Calculo do desvio total:
a1 = 0;
a2 = 0;
a3 = 0;
a4 = 0;

for s = 1:nf
    r = 4*(s-1);
    g1 = trace(G(:,:,r+1)*U1);
    g2 = trace(G(:,:,r+2)*U1);
    g3 = trace(G(:,:,r+3)*U1);
    g4 = trace(G(:,:,r+4)*U1);
    a1 = a1 + abs((p*conj(g1) + conj(p)*g1)/(2*abs(p)))^2;
    a2 = a2 + abs((p*conj(g2) + conj(p)*g2)/(2*abs(p)))^2;
    a3 = a3 + abs((p*conj(g3) + conj(p)*g3)/(2*abs(p)))^2;
    a4 = a4 + abs((p*conj(g4) + conj(p)*g4)/(2*abs(p)))^2;
    end
    %desvio durante o tempo morto:
    g3 = trace(e3);
    g4 = trace(e4);
    a5 = (p*conj(g3) + conj(p)*g3)/(2*abs(p));
    a6 = (p*conj(g4) + conj(p)*g4)/(2*abs(p));
    dev = a1 + a2 + a3 + a4 + abs(a5)^2 + abs(a6)^2;
    f0 = 1 - abs(p)/(2*S+1) + penal;
    f1 = sqrt(dev)/(2*S+1);
    F = f0 + f1; %Fidelidade
end