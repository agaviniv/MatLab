function h = blochop(E)
% Bloch sphere from quantum operation E

syms x y z
ro=0.5*eye(2)+0.5*(x*[0 1;1 0]+y*[0 -i;i 0]+z*[1 0;0 -1]);

try
    Ero = E * ro  * E'; 
catch
    disp('Dimensão da operação é diferente 2x2.')
end    

xe = trace(Ero*[0 1;1 0]);
ye = trace(Ero*[0 -i;i 0]);
ze = trace(Ero*[1 0;0 -1]);

np = 20;
[xs ys zs] = sphere(np);

xen = double(subs(xe,{x,y,z},{xs,ys,zs}));
yen = double(subs(ye,{x,y,z},{xs,ys,zs}));
zen = double(subs(ze,{x,y,z},{xs,ys,zs}));

colormap([0 0 0])
mesh(xen,yen,zen); axis([-1 1 -1 1 -1 1]); figure(gcf)
alpha(0.5); grid off; box on

