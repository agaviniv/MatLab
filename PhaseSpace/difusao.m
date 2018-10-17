function rod = difusao(acop,dq,dp,roi)
% Processo difusivo
% acop - acoplamento com o ambiente
% dq, dp - inteiros

N = size(roi,1);
I = eye(N);
U = [I(end,:);I(1:end-1,:)];
FT=sqrt(1/N)*fft(I);
V=FT*U*FT';

Dqp = U^dq * V^dp * exp(i*pi*dq*dp/N); 
rod = (1-acop)*roi + (acop/2)*(Dqp*roi*Dqp' + Dqp'*roi*Dqp);