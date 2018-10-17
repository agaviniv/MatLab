function matrix = rocoe(N,q0,p0)
% Escreve funcao de onda de estado coerente
% N e a dimensao do espaco de Hilbert

% q,p = 0,0 representa o lado esquerdo da Wigner

h=1/N;   % hbar cte de planck
dp=1/h;    % dq * dp = h/4
dq=1;    %  dq=dp (quadrado)

alpha=(q0+i*p0)/sqrt(2*h);
w0=1/N;

%%%%% gera o vetor com a funcao de onda
psi=zeros(N,1);
for j=1:N;
    q=(j-1)*dq;%-N*dq;
    p=(j-1)*dp;%-N*dp;
    psi(j)=(w0/(pi*h))^(1/4)*exp(-w0*q^2/(2*h) + sqrt(2*w0/h)*alpha*q...
                                 -alpha*conj(alpha)/2 - alpha*alpha/2);    
end

matrix = psi*(psi)' / trace(psi*(psi)') ;