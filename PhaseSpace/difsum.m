function rod = difsum(acop,dq,dp,roi)
% Processo difusivo - soma de varios deslocamentos
% acop - acoplamento com o ambiente
% dq, dp - vetore inteiros de mesmo tamanho

N = size(roi,1);
I = eye(N);
U = [I(end,:);I(1:end-1,:)];
FT=sqrt(1/N)*fft(I);
V=FT*U*FT';

M=max(size(dq)); if (M ~=max(size(dp))); error('Dimensão de dq e dp devem ser iguais'); return; end;
rod = zeros(N);
for k=1:M
    q=dq(k); p=dp(k);
    Dqp = U^q * V^p * exp(i*pi*q*p/N); 
    rod = rod + (Dqp*roi*Dqp' + Dqp'*roi*Dqp);   
end    
rod = (1-acop)*roi + (acop/(2*M))*rod;