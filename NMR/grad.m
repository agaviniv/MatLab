function ero = grad(pot,tempo,ro)

pauliZ=[1 0;0 -1];

nqbt = log(size(ro,1))/log(2); Iz=zeros(2^nqbt);
for k = 1:nqbt
    auxZ = 1;
    lb = zeros(1,nqbt);
    lb(k) = 1;
    for l = 1:nqbt
       if (lb(l) == 0)          
          auxZ = kron(auxZ,eye(2));
       elseif (lb(l) == 1)          
          auxZ = kron(auxZ,pauliZ);
       end
    end      
    Iz = Iz + auxZ;    
end

M = 100; x=1:M;
Gx = -pot*(1-2/M*x); %gradiente em forma de rampa
%Gx = pot*sin(pi/M*x);
wro = zeros(2^nqbt);
for j=1:M
    Uj = expm(i*Gx(j)*tempo*Iz);
    wro = wro + Uj*ro*Uj';
end
ero = wro/M;