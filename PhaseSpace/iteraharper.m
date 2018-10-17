function rob = iteraharper(gama,roi,nit)
% Mapa de Harper

N = size(roi,1);
I = eye(N);
U = diag(exp(-i*gama*N*cos((2*pi/N)*([0:N-1]+1/2))));
FT=sqrt(1/N)*fft(I);
V=diag(exp(-i*gama*N*cos((2*pi/N)*([0:N-1]+1/2))))*FT;

QH = U*FT'*V*FT;

for k=1:nit; rob = QH^k * roi * (QH')^k; end;