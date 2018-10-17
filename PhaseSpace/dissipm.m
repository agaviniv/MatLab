function rod = dissipm(e,a,roi)
% Modelo dissipativo
% a : <1 e >=0  

N = size(roi,1);

rod = zeros(N); G=zeros(N);
for k=0:N-1
    i = mod(-N/2+k,N); 
    ia=mod(round(i*a),N);
    P = zeros(N); 
    P(ia+1,i+1) = 1;   
    rod = rod + P*roi*P';       
    
    A = zeros(N);      I = zeros(N);
    A(ia+1,ia+1) = 1;  I(i+1,i+1) = 1;
    G=e*(G+A-I);
    
end    
rod = (1-e)*roi + e*rod; fighusi(G)