function ww = fwigsc(roin)

load mtx2qbt

wwz = zeros(9,9);

for p = 0:8

wwz(p+1,1) = (1/8)*scattcirc(roin,(V')^p*R,'z');
wwz(p+1,2) = (1/8)*scattcirc(roin,(V')^p*R*U*exp(.25*pi*i*p),'z');
wwz(p+1,3) = (1/8)*scattcirc(roin,(V')^p*R*U^2*exp(.25*pi*i*2*p),'z');
wwz(p+1,4) = (1/8)*scattcirc(roin,(V')^p*R*U^3*exp(.25*pi*i*3*p),'z');
wwz(p+1,5) = (1/8)*scattcirc(roin,(V')^p*R*U^4*exp(.25*pi*i*4*p),'z');
wwz(p+1,6) = (1/8)*scattcirc(roin,(V')^p*R*U^5*exp(.25*pi*i*5*p),'z');
wwz(p+1,7) = (1/8)*scattcirc(roin,(V')^p*R*U^6*exp(.25*pi*i*6*p),'z');
wwz(p+1,8) = (1/8)*scattcirc(roin,(V')^p*R*U^7*exp(.25*pi*i*7*p),'z');
wwz(p+1,9) = wwz(p+1,1);

end

wwy = zeros(9,9);

for p = 0:8

wwy(p+1,1) = (1/8)*scattcirc(roin,(V')^p*R,'y');
wwy(p+1,2) = (1/8)*scattcirc(roin,(V')^p*R*U*exp(.25*pi*i*p),'y');
wwy(p+1,3) = (1/8)*scattcirc(roin,(V')^p*R*U^2*exp(.25*pi*i*2*p),'y');
wwy(p+1,4) = (1/8)*scattcirc(roin,(V')^p*R*U^3*exp(.25*pi*i*3*p),'y');
wwy(p+1,5) = (1/8)*scattcirc(roin,(V')^p*R*U^4*exp(.25*pi*i*4*p),'y');
wwy(p+1,6) = (1/8)*scattcirc(roin,(V')^p*R*U^5*exp(.25*pi*i*5*p),'y');
wwy(p+1,7) = (1/8)*scattcirc(roin,(V')^p*R*U^6*exp(.25*pi*i*6*p),'y');
wwy(p+1,8) = (1/8)*scattcirc(roin,(V')^p*R*U^7*exp(.25*pi*i*7*p),'y');
wwy(p+1,9) = wwy(p+1,1);

end

www = real(wwz +wwy);


q = 0:8;  p = 0:8;
[Q P] = meshgrid(q,p);

pcolor(Q,P,www)

colormap(rwb)
shading('flat');

ww = www(1:8,1:8);