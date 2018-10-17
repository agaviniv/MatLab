function cc = concurence(roin)

yy = kron([0 -i;i 0],[0 -i;i 0]);

rotil = yy*conj(roin)*yy;

lambda = sort(eig(rotil*roin));

flambda = sqrt(lambda(4))-sqrt(lambda(3))-sqrt(lambda(2))-sqrt(lambda(1));

cc = real(max(flambda,0));