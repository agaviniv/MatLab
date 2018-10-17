%
%%

clear all

III = eye(8);
IIX = expm(-i*(pi/4)*(kron(eye(4),[0 1;1 0])));
IXI = expm(-i*(pi/4)*kron(eye(2),kron([0 1;1 0],eye(2))));
XII = expm(-i*(pi/4)*kron([0 1;1 0],eye(4)));
IIY = expm(-i*(pi/4)*kron(eye(4),[0 -i;i 0]));
IYI = expm(-i*(pi/4)*kron(eye(2),kron([0 -i;i 0],eye(2))));
YII = expm(-i*(pi/4)*kron([0 -i;i 0],eye(4)));

IXX = IXI*IIX ;
XIX = XII*IIX;
XXI = XII*IXI;
XXX = XII*IXI*IIX;

IYY = IYI*IIY;
YIY = YII*IIY;
YYI = YII*IYI;
YYY = YII*IYI*IIY;

IXY = IXI*IIY;
IYX = IYI*IIX;
XIY = XII*IIY;
YIX = YII*IIX;
XYI = XII*IYI;
YXI = YII*IXI;

XXY = XXI*IIY;
XYX = XYI*IIY;
YXX = YXI*IIX;

YYX = YYI*IIX;
YXY = YXI*IIY;
XYY = XYI*IIY;


syms x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 x13 x14 x15 x16 real
syms x17 x18 x19 x20 x21 x22 x23 x24 x25 x26 x27 x28 x29 real
syms x30 x31 x32 x33 x34 x35 x36 x37 x38 x39 x40 x41 x42 x43 real
syms x44 x45 x46 x47 x48 x49 x50 x51 x52 x53 x54 x55 x56 x57 real
syms x58 x59 x60 x61 x62 x63 x64 real
sro = [x1       x2+i*x37   x3+i*x38   x4+i*x39   x5+i*x40   x6+i*x41   x7+i*x42   x8+i*x43;
      x2-i*x37  x9         x10+i*x44  x11+i*x45  x12+i*x46  x13+i*x47  x14+i*x48  x15+i*x49;
      x3-i*x38  x10-i*x44  x16        x17+i*x50  x18+i*x51  x19+i*x52  x20+i*x53  x21+i*x54;
      x4-i*x39  x11-i*x45  x17-i*x50  x22        x23+i*x55  x24+i*x56  x25+i*x57  x26+i*x58;
      x5-i*x40  x12-i*x46  x18-i*x51  x23-i*x55  x27        x28+i*x59  x29+i*x60  x30+i*x61;
      x6-i*x41  x13-i*x47  x19-i*x52  x24-i*x56  x28-i*x59  x31        x32+i*x62  x33+i*x63;
      x7-i*x42  x14-i*x48  x20-i*x53  x25-i*x57  x29-i*x60  x32-i*x62  x34        x35+i*x64;
      x8-i*x43  x15-i*x49  x21-i*x54  x26-i*x58  x30-i*x61  x33-i*x63  x35-i*x64  x36];
    

t = 0; A=zeros(649,64);
while (t < 27)
    a = 24 * t;   
    switch (1) 
        case (t == 0)
          IN = III; 
        case (t == 1)
          IN = IIX; 
        case (t == 2)
          IN = IXI; 
        case (t == 3)
          IN = XII; 
        case (t == 4)
          IN = IIY; 
        case (t == 5)
          IN = IYI; 
        case (t == 6)
          IN = YII; 
        case (t == 7)
          IN = IXX; 
        case (t == 8)
          IN = XIX; 
        case (t == 9)
          IN = XXI; 
        case (t == 10)
          IN = XXX; 
        case (t == 11)
          IN = IYY; 
        case (t == 12)
          IN = YIY; 
        case (t == 13)
          IN = YYI; 
        case (t == 14)
          IN = YYY; 
        case (t == 15)
          IN = IXY; 
        case (t == 16)
          IN = IYX; 
        case (t == 17)
          IN = XIY;   
        case (t == 18)
          IN = YIX; 
        case (t == 19)
          IN = XYI; 
        case (t == 20)
          IN = YXI; 
        case (t == 21)
          IN = XXY; 
        case (t == 22)
          IN = XYX; 
        case (t == 23)
          IN = YXX; 
        case (t == 24)
          IN = YYX; 
        case (t == 25)
          IN = YXY; 
        case (t == 26)
          IN = XYY; 
    end

	O = IN;
	Ot = O';
    
    wro = O*sro*Ot;
    
    roaux1 = simplify(real(wro(1,5))); 
    roaux2 = simplify(imag(wro(1,5))); 
    roaux3 = simplify(real(wro(2,6))); 
    roaux4 = simplify(imag(wro(2,6))); 
    roaux5 = simplify(real(wro(3,7))); 
    roaux6 = simplify(imag(wro(3,7))); 
    roaux7 = simplify(real(wro(4,8))); 
    roaux8 = simplify(imag(wro(4,8))); 
    
    roaux9 = simplify(real(wro(1,3))); 
    roaux10 = simplify(imag(wro(1,3))); 
    roaux11 = simplify(real(wro(2,4))); 
    roaux12 = simplify(imag(wro(2,4))); 
    roaux13 = simplify(real(wro(5,7))); 
    roaux14 = simplify(imag(wro(5,7))); 
    roaux15 = simplify(real(wro(6,8))); 
    roaux16 = simplify(imag(wro(6,8)));
    
    roaux17 = simplify(real(wro(1,2))); 
    roaux18 = simplify(imag(wro(1,2))); 
    roaux19 = simplify(real(wro(3,4))); 
    roaux20 = simplify(imag(wro(3,4))); 
    roaux21 = simplify(real(wro(5,6))); 
    roaux22 = simplify(imag(wro(5,6))); 
    roaux23 = simplify(real(wro(7,8))); 
    roaux24 = simplify(imag(wro(7,8)));
        
    for j =1:64
      l = strcat('x',num2str(j));
      A(a+1,j) = double(diff(roaux1,l));
      A(a+2,j) = double(diff(roaux2,l));
	  A(a+3,j) = double(diff(roaux3,l));
      A(a+4,j) = double(diff(roaux4,l));  
      A(a+5,j) = double(diff(roaux5,l));
      A(a+6,j) = double(diff(roaux6,l));
	  A(a+7,j) = double(diff(roaux7,l));
      A(a+8,j) = double(diff(roaux8,l));  
      A(a+9,j) = double(diff(roaux9,l));
      A(a+10,j) = double(diff(roaux10,l));
	  A(a+11,j) = double(diff(roaux11,l));
      A(a+12,j) = double(diff(roaux12,l));  
      A(a+13,j) = double(diff(roaux13,l));
      A(a+14,j) = double(diff(roaux14,l));
	  A(a+15,j) = double(diff(roaux15,l));
      A(a+16,j) = double(diff(roaux16,l));
      A(a+17,j) = double(diff(roaux17,l));
      A(a+18,j) = double(diff(roaux18,l));
	  A(a+19,j) = double(diff(roaux19,l));
      A(a+20,j) = double(diff(roaux20,l));  
      A(a+21,j) = double(diff(roaux21,l));
      A(a+22,j) = double(diff(roaux22,l));
	  A(a+23,j) = double(diff(roaux23,l));
      A(a+24,j) = double(diff(roaux24,l));
    end   
   t = t + 1;
end    
A(649,1)=1; A(649,9)=1; A(649,16)=1; A(649,22)=1; A(649,27)=1; A(649,31)=1;
A(649,34)=1; A(649,36)=1;

C = zeros(64,64);
for k =1:64
    for l=1:64
        for alfa = 1:649
            C(k,l) = C(k,l)+A(alfa,k)*A(alfa,l);
        end
    end    
end

clear x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 x13 x14 x15 x16
clear x17 x18 x19 x20 x21 x22 x23 x24 x25 x26 x27 x28 x29
clear x30 x31 x32 x33 x34 x35 x36 x37 x38 x39 x40 x41 x42 x43
clear x44 x45 x46 x47 x48 x49 x50 x51 x52 x53 x54 x55 x56 x57
clear x58 x59 x60 x61 x62 x63 x64

%%
%--- Elementos de matriz referente aos espectros obtidos na tomografia ---%
rotomo = sqrt(0.5)*([1;0;0;0;0;0;0;1]*[1 0 0 0 0 0 0 1]);

t = 0;Balfa = zeros(1,649);
while (t < 27)
    a = 24 * t;   
    switch (1) 
        case (t == 0)
          IN = III; 
        case (t == 1)
          IN = IIX; 
        case (t == 2)
          IN = IXI; 
        case (t == 3)
          IN = XII; 
        case (t == 4)
          IN = IIY; 
        case (t == 5)
          IN = IYI; 
        case (t == 6)
          IN = YII; 
        case (t == 7)
          IN = IXX; 
        case (t == 8)
          IN = XIX; 
        case (t == 9)
          IN = XXI; 
        case (t == 10)
          IN = XXX; 
        case (t == 11)
          IN = IYY; 
        case (t == 12)
          IN = YIY; 
        case (t == 13)
          IN = YYI; 
        case (t == 14)
          IN = YYY; 
        case (t == 15)
          IN = IXY; 
        case (t == 16)
          IN = IYX; 
        case (t == 17)
          IN = XIY;   
        case (t == 18)
          IN = YIX; 
        case (t == 19)
          IN = XYI; 
        case (t == 20)
          IN = YXI; 
        case (t == 21)
          IN = XXY; 
        case (t == 22)
          IN = XYX; 
        case (t == 23)
          IN = YXX; 
        case (t == 24)
          IN = YYX; 
        case (t == 25)
          IN = YXY; 
        case (t == 26)
          IN = XYY; 
    end
    
    rop = IN*rotomo*IN';      % Balfa combina a ordem com roaux!!!!
    Balfa(a+1:a+24)= [real(rop(1,5)) imag(rop(1,5)) real(rop(2,6)) imag(rop(2,6))...
                      real(rop(3,7)) imag(rop(3,7)) real(rop(4,8)) imag(rop(4,8))...
                      real(rop(1,3)) imag(rop(1,3)) real(rop(2,4)) imag(rop(2,4))...
                      real(rop(5,7)) imag(rop(5,7)) real(rop(6,8)) imag(rop(6,8))...
                      real(rop(1,2)) imag(rop(1,2)) real(rop(3,4)) imag(rop(3,4))...
                      real(rop(5,6)) imag(rop(5,6)) real(rop(7,8)) imag(rop(7,8))];
    t = t + 1;                
end

Balfa(649)=0;    

%%
%--- Resolve o sistema de equacoes ---%
y = Balfa*A/C;

tomoro = [y(1)          y(2)+i*y(37)   y(3)+i*y(38)   y(4)+i*y(39)   y(5)+i*y(40)   y(6)+i*y(41)   y(7)+i*y(42)   y(8)+i*y(43);
          y(2)-i*y(37)  y(9)           y(10)+i*y(44)  y(11)+i*y(45)  y(12)+i*y(46)  y(13)+i*y(47)  y(14)+i*y(48)  y(15)+i*y(49);
          y(3)-i*y(38)  y(10)-i*y(44)  y(16)          y(17)+i*y(50)  y(18)+i*y(51)  y(19)+i*y(52)  y(20)+i*y(53)  y(21)+i*y(54);
          y(4)-i*y(39)  y(11)-i*y(45)  y(17)-i*y(50)  y(22)          y(23)+i*y(55)  y(24)+i*y(56)  y(25)+i*y(57)  y(26)+i*y(58);
          y(5)-i*y(40)  y(12)-i*y(46)  y(18)-i*y(51)  y(23)-i*y(55)  y(27)          y(28)+i*y(59)  y(29)+i*y(60)  y(30)+i*y(61);
          y(6)-i*y(41)  y(13)-i*y(47)  y(19)-i*y(52)  y(24)-i*y(56)  y(28)-i*y(59)  y(31)          y(32)+i*y(62)  y(33)+i*y(63);
          y(7)-i*y(42)  y(14)-i*y(48)  y(20)-i*y(53)  y(25)-i*y(57)  y(29)-i*y(60)  y(32)-i*y(62)  y(34)          y(35)+i*y(64);
          y(8)-i*y(43)  y(15)-i*y(49)  y(21)-i*y(54)  y(26)-i*y(58)  y(30)-i*y(61)  y(33)-i*y(63)  y(35)-i*y(64)  y(36)];


subplot(1,2,1); bar3(real(tomoro)); %axis([0.5 4.5 0.5 4.5 mn mx]);   
grid off
subplot(1,2,2); bar3(imag(tomoro)); %axis([0.5 4.5 0.5 4.5 mn mx]);
grid off      