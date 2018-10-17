function ero = atenufase(lambda,qbit,ro)
% Atenuacao de fase - Descoerencia

E{1}=[1 0;0 sqrt(1-lambda)];
E{2}=[0 0;0 sqrt(lambda)];

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

Dro=zeros(2^nqbt);
for i=1:2
        Dro = Dro + cE{i}*ro*cE{i}';
end

ero = Dro;