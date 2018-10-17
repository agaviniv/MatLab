% Tomografia de Estado para um sistema de 2 qits

%%
clear all

II = eye(4);
IX = expm(-i*(pi/4)*(kron(eye(2),[0 1;1 0])));
XI = expm(-i*(pi/4)*kron([0 1;1 0],eye(2)));
IY = expm(-i*(pi/4)*kron(eye(2),[0 -i;i 0]));
YI = expm(-i*(pi/4)*kron([0 -i;i 0],eye(2)));
XX = XI*IX;
XY = XI*IY;
YX = YI*IX;
YY = YI*IY;

syms x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 x13 x14 x15 x16 real
sro = [x1       x2+i*x11  x3+i*x12  x4+i*x13;
      x2-i*x11  x5        x6+i*x14  x7+i*x15;
      x3-i*x12  x6-i*x14  x8        x9+i*x16;
      x4-i*x13  x7-i*x15  x9-i*x16  x10];

t = 0;A=zeros(73,16);
while (t < 9)  
    a = 8 * t;   
    switch (1) 
        case (t == 0)
          IN = II; 
        case (t == 1)
          IN = IX; 
        case (t == 2)
          IN = XI; 
        case (t == 3)
          IN = XX; 
        case (t == 4)
          IN = IY; 
        case (t == 5)
          IN = YI; 
        case (t == 6)
          IN = YY; 
        case (t == 7)
          IN = XY; 
        case (t == 8)
          IN = YX; 
    end

	O = IN;
	Ot = O';
    
    wro = O*sro*Ot;
    
    roaux1 = simplify(real(wro(1,3)));  %H
    roaux2 = simplify(imag(wro(1,3))); %imaginario
    roaux3 = simplify(real(wro(2,4)));  %H
    roaux4 = simplify(imag(wro(2,4))); %imaginario
    
    roaux5 = simplify(real(wro(1,2)));  %C
    roaux6 = simplify(imag(wro(1,2))); %imaginario
    roaux7 = simplify(real(wro(3,4)));  %C
    roaux8 = simplify(imag(wro(3,4))); %imaginario
    
    for j =1:16
      l = strcat('x',num2str(j));
      A(a+1,j) = double(diff(roaux1,l));
      A(a+2,j) = double(diff(roaux2,l));
	  A(a+3,j) = double(diff(roaux3,l));
      A(a+4,j) = double(diff(roaux4,l));  
      A(a+5,j) = double(diff(roaux5,l));
      A(a+6,j) = double(diff(roaux6,l));
	  A(a+7,j) = double(diff(roaux7,l));
      A(a+8,j) = double(diff(roaux8,l));  
    end   
   t = t + 1;
end    
A(73,1)=1; A(73,5)=1; A(73,8)=1; A(73,10)=1;

C = zeros(16,16);
for k =1:16
    for l=1:16
        for alfa = 1:73
            C(k,l) = C(k,l)+A(alfa,k)*A(alfa,l);
        end
    end    
end

clear x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 x13 x14 x15 x16

%%
%--- Elementos de matriz referente aos espectros obtidos na tomografia ---%
rotomo = sqrt(0.25)*([1;1;-1;1]*[1 1 -1 1]);

t = 0;Balfa = zeros(1,73); 
while (t < 9)
    a = 8 * t;   
    switch (1) 
        case (t == 0)
          IN = II; 
        case (t == 1)
          IN = IX; 
        case (t == 2)
          IN = XI; 
        case (t == 3)
          IN = XX; 
        case (t == 4)
          IN = IY; 
        case (t == 5)
          IN = YI; 
        case (t == 6)
          IN = YY; 
        case (t == 7)
          IN = XY; 
        case (t == 8)
          IN = YX; 
    end
      
    rop = IN*rotomo*IN';      % Balfa combina a ordem com roaux!!!!
    Balfa(a+1:a+8)= [real(rop(1,3)) imag(rop(1,3)) real(rop(2,4)) imag(rop(2,4))...
                     real(rop(1,2)) imag(rop(1,2)) real(rop(3,4)) imag(rop(3,4))];
    t = t + 1;                
end

Balfa(73)=0;

%%
%--- Resolve o sistema de equacoes ---%
y = Balfa*A/C;

tomoro = [y(1)       y(2)+i*y(11)  y(3)+i*y(12)  y(4)+i*y(13);
          y(2)-i*y(11)  y(5)        y(6)+i*y(14)  y(7)+i*y(15);
          y(3)-i*y(12)  y(6)-i*y(14)  y(8)        y(9)+i*y(16);
          y(4)-i*y(13)  y(7)-i*y(15)  y(9)-i*y(16)  y(10)];

mx=max(max(real(tomoro))); mn=min(min(real(tomoro)));
mxi=max(max(imag(tomoro))); mni=min(min(imag(tomoro)));
mx=max(mx,mxi); mn=min(mn,mni);
subplot(1,2,1); bar3(real(tomoro)); axis([0.5 4.5 0.5 4.5 mn mx]);   
grid off
subplot(1,2,2); bar3(imag(tomoro)); axis([0.5 4.5 0.5 4.5 mn mx]);
grid off