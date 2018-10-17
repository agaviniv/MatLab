function fdd = fddmed(Uth,Up,nit)

dim = size(Uth,1);

A = Uth^nit;
B = Up^nit;

fdd = (norm(trace(A'*B))^2 + dim) / (dim^2 + dim);