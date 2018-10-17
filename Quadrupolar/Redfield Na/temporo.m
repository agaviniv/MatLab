function dro = temporo(ro,roeq,np,st,C,J0,J1,J2,p)
% todas as ro's devem ser matrizes densidades validas!!!
rof = ro;
dro = cell(1,np);

     %R0(1) = trace(rof - roeq);
     R0(1) = sum(diag(rof - roeq)'*[-1; 1; 1; -1]);
     R0(2) = sum(diag(rof - roeq)'*[1; 1; -1; -1]); 
     R0(3) = sum(diag(rof - roeq)'*[1; -1; 1; -1]);

for k=1:np
    dro{k} = rof;
    dt = (k)*st;
%     R0(1) = trace(rof - roeq);
%     R0(2) = sum(diag(rof - roeq)'*[-1; 1; 1; -1]);
%     R0(3) = sum(diag(rof - roeq)'*[1; 1; -1; -1]); 

    A = R0(1) * exp(-2*C*(J1+J2)*dt);
    B = R0(2) * exp(-2*C*J2*dt);
    D = R0(3) * exp(-2*C*J1*dt);
    
    % Elementos da matriz densidade
    rof(1,2) = ((1+exp(-2*C*J2*dt))*ro(1,2) + (1-exp(-2*C*J2*dt))*ro(3,4))* exp(-C*(J0+J1)*dt);
    rof(3,4) = ((1-exp(-2*C*J2*dt))*ro(1,2) + (1+exp(-2*C*J2*dt))*ro(3,4))* exp(-C*(J0+J1)*dt);
    rof(1,3) = ((1+exp(-2*C*J1*dt))*ro(1,3) + (1-exp(-2*C*J1*dt))*ro(2,4))* exp(-C*(J0+J2)*dt);
    rof(2,4) = ((1-exp(-2*C*J1*dt))*ro(1,3) + (1+exp(-2*C*J1*dt))*ro(2,4))* exp(-C*(J0+J2)*dt);
    rof(2,3) = ro(2,3) * exp(-C*(J1+J2)*dt);
    rof(1,4) = ro(1,4) * exp(-C*(J1+J2)*dt);
    rof(1,1) = p -(1/4)*(A -B -D);
    rof(2,2) = 2*p/3 +(1/4)*(+A +B -D);
    rof(3,3) = p/3 +(1/4)*(+A -B +D);
    rof(4,4) = 0*p -(1/4)*(A +B +D);
            rof(2,1)=rof(1,2); rof(3,1)=rof(1,3); rof(4,1)=rof(1,4);
            rof(3,2)=rof(2,3); rof(4,2)=rof(2,4); rof(4,3)=rof(3,4);
    
   % dro{k} = rof;
end