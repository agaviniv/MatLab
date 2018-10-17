function ero = acoplaJ(J,tempo,qbits,ro)

pauliZ=[1 0;0 -1];

tabq = ['a' 'b' 'c' 'd' 'e' 'f' 'g' 'h' 'i' 'j' 'k' 'l']; % max 12 qbits

nqbt = log(size(ro,1))/log(2); Iz=eye(2^nqbt);
s = size(qbits,2);
for j = 1:s    
    auxZ = 1;
    lb = (tabq == qbits(j));
    for l = 1:nqbt
       if (lb(l) == 0)          
          auxZ = kron(auxZ,eye(2));
       elseif (lb(l) == 1)          
          auxZ = kron(auxZ,pauliZ);
       end
    end      
    Iz = Iz * 0.5*auxZ;        
end

Uj = expm(i*2*pi*(J*tempo)*Iz);
ero = Uj*ro*Uj';