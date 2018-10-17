function fun = lorentzd(w,w0,wnut,m0,T1,T2)
%Dispersao positiva na direcao -x

fun = m0*(wnut*T2^2*(w-w0))./(1 + T2^2*(w-w0).^2 + wnut^2*T1*T2);