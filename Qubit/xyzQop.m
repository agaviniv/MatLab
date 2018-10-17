% Bloch sphere from quantum operation

clear all

syms x y z
ro=0.5*eye(2)+0.5*(x*[0 1;1 0]+y*[0 -i;i 0]+z*[1 0;0 -1]);

E = sqrt(0.5)*[1 0;1 2*sqrt(0.5)];
%Ero = E * ro  * E'; 
Ero = atenuamp(1,.5,'a',ro);

xe = trace(Ero*[0 1;1 0]);
ye = trace(Ero*[0 -i;i 0]);
ze = trace(Ero*[1 0;0 -1]);

np = 15;
[xs ys zs] = sphere(np);

xen = double(subs(xe,{x,y,z},{xs,ys,zs}));
yen = double(subs(ye,{x,y,z},{xs,ys,zs}));
zen = double(subs(ze,{x,y,z},{xs,ys,zs}));

colormap([0 0 0])
subplot(1,2,1); mesh(xs,ys,zs); axis([-1 1 -1 1 -1 1]); alpha(0.5); grid off
subplot(1,2,2); mesh(xen,yen,zen); axis([-1 1 -1 1 -1 1]); alpha(0.5); grid off



roI = 0.5*eye(2);
for kk=1:40
    SroI = atenuamp(1,.05,'a',roI);
    eta(kk) = 2*trace((SroI - roI)^2);
    roI = SroI;
end

figure;
plot(eta)
