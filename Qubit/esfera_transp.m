function valor = esfera_transp(raio, npontos);
%
%
%  Desenha uma esfera transparente.
%
%
%

n_theta = 0.1;

n_phi = 0.25;

hold on;

for ttt = 0:0.25:1
   
   mm = 1;
   
   for ppp = 0:1/200:2
      
      theta = ttt * pi;  phi = ppp * pi;
      
      xx(mm) = raio * cos(phi) * sin(theta);
      
      yy(mm) = raio * sin(phi) * sin(theta);

		zz(mm) = raio * cos(theta);
      
      mm = mm + 1;   
      
   end;
   
   plot3(xx,yy,zz,'-k');

   
end;


for ppp = 0:0.1:2
   
   mm = 1;
   
   for ttt = 0:1/200:1
      
      theta = ttt * pi;  phi = ppp * pi;
      
      xx(mm) = raio * cos(phi) * sin(theta);
      
      yy(mm) = raio * sin(phi) * sin(theta);

		zz(mm) = raio * cos(theta);
      
      mm = mm + 1;   
      
   end;
   
   plot3(xx,yy,zz,'-k');

   
end;


hold off;

