function rob = iterabaker(roi,nit)
% Mapa do padeiro

N = size(roi,1);
I = eye(N);
U = [I(end,:);I(1:end-1,:)];
FT=sqrt(1/N)*fft(I);
R=FT^2;
V=FT*U*FT';
FTN2=sqrt(2/N)*kron(eye(2),fft(eye(N/2)));

QB = inv(FT)*FTN2;

for k=1:nit; rob = QB^k * roi * (QB')^k; end;