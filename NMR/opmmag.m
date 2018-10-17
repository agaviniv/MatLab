function I_c = opmmagI(s,c);
%
%  Calcula as matrizes Ix, Iy e Iz (em unidades de h / 2*pi)  
%
%  s e o valor do Spin (Ex: 1/2, 3/2, 1, ...)
%  c e a componente do operador de spin 

mI = s:-1:-s;  
syms a;

%---------------------- Calculatsg the marix Ix, Iy and Iz ------------------------------

for ls=1:1:length(mI)
   for col=1:1:length(mI)
      
      if ls == col dmm = 1; else dmm = 0; end;
      
      if ls == col-1 dmp = 1; else dmp = 0; end;
      
      if ls == col+1 dpm = 1; else dpm = 0; end;
      
  		Ip = sqrt((s - mI(col)) * (s + mI(col) + 1)) * dmp;
      
      Im = sqrt((s + mI(col)) * (s - mI(col) + 1)) * dpm;
          
      if (c == 'x')
        I_c(ls, col) = (Ip + Im)/2; 
      elseif (c == 'y')
        I_c(ls, col) = (Ip - Im)/(2*i); 
      elseif (c == 'z')
        I_c(ls, col) = mI(col) * dmm;  
      end;  

   end;
end;