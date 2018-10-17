
clear all

psi = sqrt(0.5)*[1;1]; %psi=[0;1]
%|ABC> : A - amplitude, B - fase e C - sistema
ro = kron([1 0;0 0],kron([1 0;0 0],psi*psi')); 

xac = kron([1 0;0 0],eye(4))+kron([0 0;0 1],kron(eye(2),[0 1;1 0]));

np = 200;
at = 10000e-3; t1=1; t2=0.05;

Iz=zeros(1,np);Ix=zeros(1,np);dd=zeros(np,8);cc=zeros(np,6);
mm=zeros(np,3); ros = cell(1,np);
for k = 0:np-1
    t = k*at/np;  
    lambda = (asin((1-exp(-t/t2))^0.5));
    theta = (asin((1-exp(-t/t1))^0.5));
    
    %Rotacao sobre um qbit 
    ry = @(th) expm(-i*th/2*[0 -i;i 0]);

    cry_b =  kron(eye(2^2),[1 0;0 0]) ...
           + kron(eye(2),kron(ry(lambda),[0 0;0 1]));

    cry_a =  kron(eye(2^2),[1 0;0 0]) ...
           + kron(ry(theta),kron(eye(2),[0 0;0 1]));
    %%%%%%%%%%% 
             rr{k+1} = xac*cry_a*cry_b * ro * (xac*cry_a*cry_b)';   
    %%%%%%%%%%%
    ros = rr{k+1};    dd(k+1,:) = diag(ros);
    cc(k+1,:) = [ros(1,2) ros(3,4) ros(1,3) ros(2,4) ros(2,3) ros(1,4)];
    
    roa = trmsq(trmsq(ros)); mm(k+1,:) = [roa(1,1) roa(2,2) roa(1,2)];
    
    Ix(k+1) = trace(roa*[0 1;1 0]);      % FID
    Iz(k+1) = trace(roa*[1 0;0 -1]);     % Magnetizacao longitudinal
    
end

plot((0:np-1)*at/np,[Ix;Iz])
