% Program to homonuclear system
%  two spin 1/2
%
clear all

%Pauli's matrices
x1=kron([0 1;1 0],eye(2));
x2=kron(eye(2),[0 1;1 0]);

y1=kron([0 -i;i 0],eye(2));
y2=kron(eye(2),[0 -i;i 0]);

z1=kron([1 0;0 -1],eye(2));
z2=kron(eye(2),[1 0;0 -1]);


%#########################################################################
%homonuclear system
roeq=z1+z2;
frq=200;            %dibromothiophene cs=161.63, J=5.77 PRA(67),042322
J=100;
h = pi*(frq*z1 + 0.5*J*z1*z2);
T2=0.6;


%#########################################################################
%pulse timing
roin=(y1+y2); %pi/2 x:pulse in roeq


%#########################################################################
%acquisition time
at=3;  %seconds
np=8800;

dt=at/np;
Io = y1+y2;
fidt=zeros(1,np);
for k=1:np
    U= expm(-i* h * (k*dt));
    rot = U*roin*U';
   % rot=atenufase(exp(-k*dt/T2),'a',rot);rot=atenufase(exp(-k*dt/T2),'b',rot);
    fidt(k) = trace(rot*Io)*exp(-k*dt/T2);
end    


%#########################################################################
%processing
timeax = dt*[1:np];
plot(timeax,real(fidt),'', 'DisplayName', 'fidt', 'YDataSource', 'fidt'); figure(gcf)


nps=1000; 
sw = 1/dt;
freqax = sw*[-nps/2:(nps-1)/2]/nps;
spc1 = fft(fidt,nps);  spc2 = fftshift(spc1);


figure;
plot (freqax,real(spc2), 'DisplayName', 'spc2', 'YDataSource', 'spc2'); figure(gcf)