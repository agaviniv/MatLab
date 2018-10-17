
clear all

ro = 0.5*diag([1 0 0 1]);
%ro = (1/2)*[1;0;0;1]*[1 0 0 1];


%Wigner em funcao de p
W0p = @(p) (1/8)*(ro(1,1)+(-1)^p*ro(3,3)+2*real(ro(2,4)*exp(i*pi*p/2)));
W1p = @(p) (1/4)*(real(ro(1,2)*exp(-i*pi*p/4))+real(ro(3,4)*exp(i*3*pi*p/4)));
W2p = @(p) (1/8)*(ro(2,2)+(-1)^p*ro(4,4)+2*real(ro(1,3)*exp(-i*pi*p/2)));
W3p = @(p) (1/4)*(real(ro(2,3)*exp(-i*pi*p/4))+real(ro(1,4)*exp(-i*3*pi*p/4)));

W_p = [W0p(3) W1p(3) W2p(3) W3p(3);
       W0p(2) W1p(2) W2p(2) W3p(2);
       W0p(1) W1p(1) W2p(1) W3p(1);
       W0p(0) W1p(0) W2p(0) W3p(0)]; 

%Wigner em funcao de q
Wq0 = @(q) (1/8)*(ro(mod(q,4)+1,1)+ro(mod(q-1,4)+1,2)+ro(mod(q-2,4)+1,3)+ro(mod(q-3,4)+1,4));
Wq1 = @(q) (exp(i*pi*q/4)/8)*(ro(mod(q,4)+1,1)-ro(mod(q-2,4)+1,3)+i*(ro(mod(q-3,4)+1,4))-ro(mod(q-1,4)+1,2));
Wq2 = @(q) (exp(i*pi*q/2)/8)*(ro(mod(q,4)+1,1)-ro(mod(q-1,4)+1,2)+ro(mod(q-2,4)+1,3)-ro(mod(q-3,4)+1,4));
Wq3 = @(q) (exp(i*3*pi*q/4)/8)*(ro(mod(q,4)+1,1)-ro(mod(q-2,4)+1,3)-i*(ro(mod(q-3,4)+1,4))-ro(mod(q-1,4)+1,2));

W_q = [Wq3(0) Wq3(1) Wq3(2) Wq3(3);
       Wq2(0) Wq2(1) Wq2(2) Wq2(3);
       Wq1(0) Wq1(1) Wq1(2) Wq1(3);
       Wq0(0) Wq0(1) Wq0(2) Wq0(3)];