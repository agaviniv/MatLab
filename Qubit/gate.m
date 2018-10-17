function ero = gate(porta,qbit,ro)
% ''porta'' sobre um qbit
   
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
      aux = kron(aux,porta);      
   end
end      

U = aux;

try
    ero = U*ro*U';
catch
    ero = nan;
    disp('Dimensão da porta não compatível com a dimensão de ro.')
end    