function mm = scattcirc(roin,opu,sigma)

dim = size(roin,2);

ha = sqrt(.5)*kron([1 1;1 -1],eye(dim));

ctru = eye(2*dim);

ll = zeros(1,2*dim);

ctru(dim+1,dim+1:end) = opu(1,:);
ctru(dim+2,dim+1:end) = opu(2,:);
ctru(dim+3,dim+1:end) = opu(3,:);
ctru(dim+4,dim+1:end) = opu(4,:);

ro = kron(diag([1 0]),roin);

U = ha*ctru*ha;

ro = U*ro*U';

if sigma == 'x'
    m = kron([0 1;1 0],eye(dim));
elseif sigma == 'y'
    m = kron([0 -i;i 0],eye(dim));
elseif sigma == 'z'
    m = kron([1 0;0 -1],eye(dim));
end   

mm = trace(ro*m);