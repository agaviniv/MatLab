function matW = calcwig(roh);
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
mm= max(max(W)); W(p+1,:) = -mm; W(:,q+1) = -mm;
%W(p+1,:) = -W(p+1,:); W(:,q+1) = -W(:,q+1);
%W(p+1,:) = -1/(N^2); W(:,q+1) = 1/(N^2); % Normalizacao
%W(p+1,:) = W(1,:); W(:,q+1) = W(:,1);

