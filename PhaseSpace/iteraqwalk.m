function roqw = iteraqwalk(roi,nit)
% Mapa do "quantum walk"

N = size(roi,1)/2;
I = eye(N);
U = kron([I(end,:);I(1:end-1,:)]^(-1),diag([1 0]))+...
    kron([I(end,:);I(1:end-1,:)],diag([0 1]));

QW = U * kron(eye(N),sqrt(1/2)*[1 1;1 -1]);

for k=1:nit; ro = QW^k * roi * (QW')^k; end;

roqw = trlsq(ro);