% Defini matrizes para estudo da representacao do estado
%           quantico no espaco de fases discreto

N = input('Dimensão(2^n):');
I = eye(N);
U = [I(end,:);I(1:end-1,:)];
FT=sqrt(1/N)*fft(I);
R=FT^2;
V=FT*U*FT';
FTN2=sqrt(2/N)*kron(eye(2),fft(eye(N/2)));

T=cell(N,N); A=cell(N,N);
for q=0:N-1
    for p=0:N-1
        T{q+1,p+1}=U^q*V^p*exp(i*pi*q*p/N);
        A{q+1,p+1}=(0.5/N)*U^q*R*V^(-p)*exp(i*pi*q*p/N);
    end    
end    
clear q p
