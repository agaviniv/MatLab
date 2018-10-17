function mm = mprobqp(roh);
% Programa para o calculo da Funcao de Wigner

N = size(roh,1);   % Dimensao do Espaco de Hilbert do Sistema Alvo: 2^(nº de qubits)

S = 0;
for q = 0:1:(2*N)
   for p = 0:1:(2*N)
      for n = 0:1:(N-1)
         a = 1+mod(q-n,N);
         b = 1+mod(n,N);
         S = S + roh(a,b) * exp((i)*(pi/N)*p*(q - 2*n));
      end
      W(p+1,q+1) = real(S /(2*N));
      S = 0;          
   end
end   
matW = W(1:end-1,1:end-1);
subplot(1,2,1); bar(sum(matW)); axis([0 2*N 0 1]); xlabel('q')
subplot(1,2,2); bar(sum(matW')); axis([0 2*N 0 1]); xlabel('p')
