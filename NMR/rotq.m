function ero = rotq(flip,fase,qbit,ro)
% Rotacao - Pulso "hard"

pauli{1}=[0 1;1 0];
pauli{2}=[0 -i;i 0];
pauli{3}=[1 0;0 -1];

nqbt = log(size(ro,1))/log(2);
tabq = ['a' 'b' 'c' 'd' 'e' 'f' 'g' 'h' 'i' 'j' 'k' 'l']; % max 12 qbits

aux1 = 1;
aux2 = 1;
aux3 = 1;
I = cell(1,3); I{1}=zeros(2^nqbt);I{2}=zeros(2^nqbt);I{3}=zeros(2^nqbt);

s = size(qbit,2);
for k = 1:s
    aux1 = 1;
    aux2 = 1;
    aux3 = 1;
    lb = (tabq == qbit(k));
    for l = 1:nqbt
       if (lb(l) == 0)
          aux1 = kron(aux1,eye(2));
          aux2 = kron(aux2,eye(2));
          aux3 = kron(aux3,eye(2));
       elseif (lb(l) == 1)
          aux1 = kron(aux1,pauli{1});
          aux2 = kron(aux2,pauli{2});
          aux3 = kron(aux3,pauli{3});
       end
    end      
    I{1} = I{1} + aux1;
    I{2} = I{2} + aux2;
    I{3} = I{3} + aux3;    
end

Ip = cos(fase)*I{1}-sin(fase)*I{2};
U = expm(-i*0.5*flip*Ip);

ero = U*ro*U';