function Zout = ncalcZ(Uin)

N = size(Uin,1);

I = eye(N);
U = [I(end,:);I(1:end-1,:)];
FT=sqrt(1/N)*fft(I);
R=FT^2;
V=FT*U*FT';

A=cell(N,N);
for q=0:N-1
    for p=0:N-1
        A{q+1,p+1}=(0.5/N)*U^q*R*V^(-p)*exp(i*pi*q*p/N);
    end    
end 

k=1; Z = zeros(1,4*N^2);
for q=0:N-1
    for p=0:N-1
        for nq=0:N-1
            for np=0:N-1 
                Z(k)=N*real(trace(A{q+1,p+1}*Uin*A{nq+1,np+1}*Uin')); 
                k = k + 1;
            end
        end
    end
end    
Zout = reshape(Z,sqrt(size(Z,2)),sqrt(size(Z,2))); %Zout(p,q) -> p:lin,q:col