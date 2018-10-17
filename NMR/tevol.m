function ket = tevol(ketin,H,dt,tevol)
%tevol(ketin,H,dt,tevol) -> Evolui no ketin sobre a acao de H ate o tempo
%tevol; dt e o incremento temporal

sx = [0 1;1 0]; sy = [0 -i;i 0]; sz = [1 0;0 -1];

np = round(tevol/dt);

vb = zeros(np,3);
vb(1,:) = [ketin'*sx*ketin ketin'*sy*ketin ketin'*sz*ketin];
for kk=2:np
    U = expm(-i*H*dt*kk);
    ket = U*ketin;
    vb(kk,:) = [ket'*sx*ket ket'*sy*ket ket'*sz*ket];
end    

%[X,Y,Z] = sphere(50); surf(.999*X,.999*Y,.999*Z); 

colormap([.5 .5 0.2]);
hold on
p=plot3(vb(:,1),vb(:,2),vb(:,3),'.');
esfera_transp(.999,40);
gp=get(p); s=set(p,'LineWidth',2);
view(3)