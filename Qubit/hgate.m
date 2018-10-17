function ero = hgate(qbit,ro)
% Operacao de Hadamard sobre um qbit

h = sqrt(0.5)*[1 1;1 -1];

nqbt = log(size(ro,1))/log(2);
tabq = ['a' 'b' 'c' 'd' 'e' 'f' 'g' 'h' 'i' 'j' 'k' 'l']; % max 12 qbits

aux = 1;

lb = (tabq==qbit);
if (lb==zeros(1,size(tabq,1)))
    disp('O qbit escolhido não está definido.')
    return;
end 
for l = 1:nqbt
   if (lb(l) == 0)
      aux = kron(aux,eye(2));      
   elseif (lb(l) == 1)
      aux = kron(aux,h);      
   end
end      

H = aux;

ero = H*ro*H';