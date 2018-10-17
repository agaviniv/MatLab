function ero = bitflip(p,qbit,ro)
% Atenuacao de amplitude - Dissipacao de energia

E{1}=sqrt(p)*eye(2);
E{2}=sqrt(1-p)*[0 1;1 0];

nqbt = log(size(ro,1))/log(2);
tabq = ['a' 'b' 'c' 'd' 'e' 'f' 'g' 'h' 'i' 'j' 'k' 'l']; % max 12 qbits

aux0 = 1;
aux1 = 1;
lb = (tabq==qbit);
for l = 1:nqbt
   if (lb(l) == 0)
      aux0 = kron(aux0,eye(2));
      aux1 = kron(aux1,eye(2));     
   elseif (lb(l) == 1)
      aux0 = kron(aux0,E{1});
      aux1 = kron(aux1,E{2});     
   end
end      
cE{1} = aux0;
cE{2} = aux1;
%Ordem das celulas: quanto menor o indice -> qbit menos significativo

sro=zeros(2^nqbt);
for m=1:2
        sro = sro + cE{m}*ro*cE{m}';
end

ero = sro;