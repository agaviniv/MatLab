function wro = trlsq(roin)
% Traco parcial sobre o qbit menos significante
%       (Least Significant Qubit)
% Ex. Trlsq(ro_AB)= rot_A

N = size(roin,1);

sp0 = kron(eye(N/2),[1 0]);
sp1 = kron(eye(N/2),[0 1]);

wro = sp0*roin*sp0' + sp1*roin*sp1'; 