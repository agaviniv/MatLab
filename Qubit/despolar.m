function ero = despolar(p,qbit,ro)
% Atenuacao de amplitude - Dissipacao de energia

E{1}=sqrt(1-p)*eye(2);
E{2}=sqrt(p/3)*[0 1;1 0];
E{3}=sqrt(p/3)*[0 -i;i 0];
E{4}=sqrt(p/3)*[1 0;0 -1];

nqbt = log(size(ro,1))/log(2);
tabq = ['a' 'b' 'c' 'd' 'e' 'f' 'g' 'h' 'i' 'j' 'k' 'l']; % max 12 qbits

aux0 = 1;
aux1 = 1;
aux2 = 1;
aux3 = 1;
lb = (tabq==qbit);
for l = 1:nqbt
   if (lb(l) == 0)
      aux0 = kron(aux0,eye(2));
      aux1 = kron(aux1,eye(2));
      aux2 = kron(aux2,eye(2));
      aux3 = kron(aux3,eye(2));
   elseif (lb(l) == 1)
      aux0 = kron(aux0,E{1});
      aux1 = kron(aux1,E{2});
      aux2 = kron(aux2,E{3});
      aux3 = kron(aux3,E{4});
   end
end      
cE{1} = aux0;
cE{2} = aux1;
cE{3} = aux2;
cE{4} = aux3;%Ordem das celulas: quanto menor o indice -> qbit menos significativo

sro=zeros(2^nqbt);
for m=1:4
        sro = sro + cE{m}*ro*cE{m}';
end

ero = sro;