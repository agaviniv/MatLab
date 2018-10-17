function fun = lorentza(w,w0,wnut,m0,T1,T2)
%Absorcao positiva na direcao y

fun = -m0*(wnut*T2)./(1 + T2^2*(w-w0).^2 + wnut^2*T1*T2);