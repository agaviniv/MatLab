function ero = tgate(qbit,ro)
% Operacao de Hadamard sobre um qbit

t = [1 0;0 exp(i*pi/4)];

nqbt = log(size(ro,1))/log(2);
tabq = ['a' 'b' 'c' 'd' 'e' 'f' 'g' 'h' 'i' 'j' 'k' 'l']; % max 12 qbits

aux = 1;

lb = (tabq==qbit);
for l = 1:nqbt
   if (lb(l) == 0)
      aux = kron(aux,eye(2));      
   elseif (lb(l) == 1)
      aux = kron(aux,t);      
   end
end      

T = aux;

ero = T*ro*T';