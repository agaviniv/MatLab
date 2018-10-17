function [Ix_ss Iy_ss Iz_ss] = mangqbts(nqbt)

pauli{1} = 2*opmmag(1/2,'x');
pauli{2} = 2*opmmag(1/2,'y');
pauli{3} = 2*opmmag(1/2,'z');

for k = 1:nqbt
    auxX = 1;
    auxY = 1;
    auxZ = 1;
    lb = zeros(1,nqbt);
    lb(k) = 1;
    for l = 1:nqbt
       if (lb(l) == 0)
          auxX = kron(auxX,eye(2));
          auxY = kron(auxY,eye(2));
          auxZ = kron(auxZ,eye(2));
       elseif (lb(l) == 1)
          auxX = kron(auxX,pauli{1});
          auxY = kron(auxY,pauli{2});
          auxZ = kron(auxZ,pauli{3});
       end
    end      
    Ix_ss{k} = 0.5*auxX;
    Iy_ss{k} = 0.5*auxY;
    Iz_ss{k} = 0.5*auxZ;
end
