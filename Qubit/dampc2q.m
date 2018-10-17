function ro=dampc2q(sa,sb,ea,eb,ela,elb)

psi=kron(sa,kron(sb,kron(ea,kron(eb,kron(ela,elb)))));


%parametros de atenuacao [pda pdb ada adb]
parat = [0.08 0.08 0.8 0.8];

np = 40;
dt = 1e-4
for k = 0:np-1
    t = k*dt;  
    theta = (asin((1-exp(-t/p))^2));
    %Rotacao controlada da atenuacao de fase 
    ry = @(theta) expm(-i*theta*[0 -i;i 0]);

    cry_ea = kron([1 0;0 0],eye(2^5)) ...
             + kron([0 0;0 1],kron(eye(2),kron(ry(parat(1)),eye(2^3)));

    cry_eb = kron(eye(2),kron([1 0;0 0],eye(2^4))) ...
             + kron(eye(2),kron([0 0;0 1],kron(eye(2),kron(ry(parat(2)),eye(2^2)))));

    psit = cry_eb * cry_ea * psi;
    
    %Medidas projetivas - estados reduzidos
    prj0_ea = kron(eye(2^2),kron([1 0],eye(2^3)));
    prj1_ea = kron(eye(2^2),kron([0 1],eye(2^3)));
    
    psi_reda = prj0_ea*psi + prj1_ea*psi;
    
    prj0_eb = kron(eye(2^3),kron([1 0],eye(2^2)));
    prj1_eb = kron(eye(2^3),kron([0 1],eye(2^2)));
    
    psi_redb = prj0_eb*psi + prj1_eb*psi;
    
    %qual o estado de ea e eb????

end                