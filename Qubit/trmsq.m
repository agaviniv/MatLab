function wro = trmsq(roin)
% Traco parcial sobre o qbit mais significante
%       (Most Significant Qubit)
% Ex. Trmsq(ro_AB)= rot_B

N = size(roin,1);

sp0 = kron([1 0],eye(N/2));
sp1 = kron([0 1],eye(N/2));

wro = sp0*roin*sp0' + sp1*roin*sp1'; 

