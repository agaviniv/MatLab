function matrix = fighusi(ro)
%%%%%%%%%%% Escreve funcao de onda de estado coerente

N=size(ro,1);
h=2.0^5;      %hbar cte de planck


%sumpsi=zeros(N,1);
for q0 = 0:N
    for p0 = 0:N
        dp=1;
        dq=1;

        alpha=(q0+i*p0)/sqrt(2*h);
        w0=1.0/32;
        wpsi=zeros(N,1);

        %%%%% gera o vetor com a funcao de onda
        for j=1:N;
            q=(j-1)*dq-N*dq/2;
            p=(j-1)*dp-N*dp/2;
            wpsi(j)=(w0/(pi*h))^(1/4)*exp(-w0*q^2/(2*h)+sqrt(2*w0/h)*alpha*q-alpha*conj(alpha)/2-alpha*alpha/2);    
        end
%         sumpsi = sumpsi + wpsi;
        
        husimi(p0+1,q0+1)=wpsi'*ro*wpsi;
        
    end
end    
husimi(N+1,N+1)=husimi(1,1);

%Husimiop = sumpsi*sumpsi';

xe = [0:1:N]; ye = [0:1:N];         % os eixos 0, 1, 2, 3, ..., 8 (paradois q-bits)

pcolor(xe, ye, real(husimi));      % desenha os quadrados

colormap(summer);                 % esquema de cores

shading('flat');    			% retira as linhas da grade

