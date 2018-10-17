function matriz = rocoerente(N,q0,p0)
%%%%%%%%%%% Escreve funcao de onda de estado coerente
% escrita por J. Paulo / CBPF

% N=2^6;
% q0=0;
% p0=50;

h=1;      %hbar cte de planck
alpha=(q0+i*p0)/sqrt(2*h);
w0=1;
dp=1;
dq=1;

psi=zeros(N,1);

%%%%% gera o vetor com a funcao de onda
for j=1:N;
    q=(j-1)*dq-N*dq/2;
    p=(j-1)*dp-N*dp/2;
    psi(j)=(w0/(pi*h))^(1/4)*exp(-w0*q^2/(2*h)+sqrt(2*w0/h)*alpha*q...
                                 -alpha*conj(alpha)/2-alpha*alpha/2);                             
end
n = sqrt(psi'*psi); psi=psi/n;
matriz = psi*psi';
