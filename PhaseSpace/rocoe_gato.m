%%%%%%%%%%% Escreve funcao de onda de estado coerente

N=2^6;
h=0.25;      %hbar cte de planck

q0=4;
p0=16;
dp=1;
dq=1;

alpha=(q0+i*p0)/sqrt(2*h);
w0=1.0;
psi=zeros(N,1);

%%%%% gera o vetor com a funcao de onda
for j=1:N;
    q=(j-1)*dq-N*dq/2;
    p=(j-1)*dp-N*dp/2;
    psi(j)=(w0/(pi*h))^(1/4)*exp(-w0*q^2/(2*h)+sqrt(2*w0/h)*alpha*q-alpha*conj(alpha)/2-alpha*alpha/2);    
end


q0=16;
p0=16;
dp=1;
dq=1;

alpha=(q0+i*p0)/sqrt(2*h);
w0=1;
psi1=zeros(N,1);

%%%%% gera o vetor com a funcao de onda
for j=1:N;
    q=(j-1)*dq-N*dq/2;
    p=(j-1)*dp-N*dp/2;
    psi1(j)=(w0/(pi*h))^(1/4)*exp(-w0*q^2/(2*h)+sqrt(2*w0/h)*alpha*q-alpha*conj(alpha)/2-alpha*alpha/2);    
end

psig = psi + psi1;
figwig(psig*psig')