function wro = wig2ro(fwig)

N = 0.5*size(fwig,1);

I = eye(N);
U = [I(end,:);I(1:end-1,:)];
FT=sqrt(1/N)*fft(I);
R=FT^2;
V=FT*U*FT'; %definindo V=FT'*U*FT mantem-se o expoente negativo de V em A

wro=zeros(N);
for q=0:2*N-1
     for p=0:2*N-1         
         A=(0.5/N)*U^q*R*V^(p)*exp(i*pi*q*p/N); %sinal do expoente de V
         wro = wro+ N*fwig(p+1,q+1)*A;         %trocado em relação a psop  
     end    
end
