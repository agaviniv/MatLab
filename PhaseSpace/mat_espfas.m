% Calcula as Matrizes do Espaço de Fases
%

%clear all
 
D = input('dim:');%3; % dim - 1
N = D - 1;

U = zeros(D,D); V = zeros(D,D); R = zeros(D,D);
FT = zeros(D,D); A = zeros(D,D); FTN_2 = zeros(D,D);
QB = zeros(D,D); opr = zeros(D,D);

% Operadores de translacao e reflexao
for m = 0:1:N
  for n = 0:1:N 
    if (m == mod(n+1,N+1))
       U(m+1,n+1) = 1;
    else
       U(m+1,n+1) = 0;
    end
    if (m == mod(N+1-n,N+1))
       R(m+1,n+1) = 1;
    else
       R(m+1,n+1) = 0;
    end
    if (m == n)
       V(m+1,n+1) = exp(-i*2*pi*n/(N+1));
    else
       V(m+1,n+1) = 0;
    end

  end 
end


% Transformada de Fourier
for l = 0:N
    FT(:,l+1) = sqrt(1/D)*(diag(V)).^(l);
end    


% Operadores de pontos do espaco de fases
l=1;
for q = 0:1:N
   for p = 0:1:N
      opr = exp(i*(pi)*q*p/(N+1)) * U^q * R * V^(-p);
      A(1:N+1,1:N+1,l)= opr/(N+1);
      l = l+1;
   end
end   


%/////////\\\\\\\\\
% Mapa do Padeiro
FTN_2 = sqrt(2/D)*kron(eye(2),fft(eye(D/2)));
QB = inv(FT)*FTN_2;

