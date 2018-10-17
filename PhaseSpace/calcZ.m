function Zout = calcZ(Uin)

N = size(Uin,1);

I = eye(N);
U = [I(end,:);I(1:end-1,:)];
FT=sqrt(1/N)*fft(I);
R=FT^2;
V=FT*U*FT';

k=1; Z = zeros(1,4*N^2);
for q=0:N-1
    for p=0:N-1
        Aqp=(0.5/N)*U^q*R*V^(-p)*exp(i*pi*q*p/N);
        for nq=0:N-1
            for np=0:N-1
                Anqp=(0.5/N)*U^nq*R*V^(-np)*exp(i*pi*nq*np/N);
                Z(k)=N*real(trace(Aqp*Uin*Anqp*Uin')); 
                k = k + 1;
            end
        end
    end
end    
Zout = reshape(Z,sqrt(size(Z,2)),sqrt(size(Z,2))); %Zout(p,q) -> p:lin,q:col