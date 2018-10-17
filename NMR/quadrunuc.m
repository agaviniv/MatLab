clear all

%Pauli's matrices
x=opmmag(3/2,'x');
y=opmmag(3/2,'y');
z=opmmag(3/2,'z');


%#########################################################################
%quadrupolar system set-up
N = 4;
alfa = 4e-8;
frqp=1.6e3; frql = 100e6;
h = -2*pi*frql*z + (pi/3)*frqp*(3*z^2 - 3/2*(3/2+1)*eye(4));
roeq=(1/N)*eye(4)-alfa*h; roeq = diag([1/2 1/3 1/6 0]);
T2 = 0.5;
T1a = 1.8; T1b = 1.8;



%       @@@@@@@@@    EXPERIMENT' SEQUENCE SIMULATION   @@@@@@@@@
%#########################################################################
%pulse initialization - SMP (ideal)
P1 = [1 0 0 0;0 0 0 1;0 1 0 0;0 0 1 0];
P2 = [1 0 0 0;0 0 1 0;0 0 0 1;0 1 0 0];
roin=(roeq+P1*roeq*P1'+P2*roeq*P2')/3;

ro = cell(1,4);
for kk = 0:40
%#########################################################################
%free evolution time
dt = kk*0.001;

U0 = expm(-i*h*dt);
ro0 = U0*roin*U0';

lamba = 0.05;%1 - exp(-dt/T2); 
lambb = 0.05;%1 - exp(-dt/T2);
ro{1} = atenufase(lamba,'a',ro0);
ro{2} = atenufase(lambb,'b',ro{1});

gamaa = 1 - exp(-dt/T1a); gamab = 1 - exp(-dt/T1b);
ro{3} = atenuamp(2/3,gamaa,'a',ro{2});
ro{4} = atenuamp(5/6,gamab,'b',ro{3});

dd(kk+1,:) = real(diag(ro{4}));
roin=ro{4};
%#########################################################################
%tomographic processing

end