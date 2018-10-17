function estimated_gama=prefit(np,p,start,curva,curvafun)

% Fitting com canal de atenuacao de amplitude
t = 0:np-1; 
figure;
plot(t,curva,'ro'); hold on; h = plot(t,curva,'b'); hold off;
title('Input data'); ylim([0 0.7])


options = optimset('Simplex','on','TolX',0.001);
estimated_gama = fminsearch(@(x)disalfafit(x,t,curva,p,h),start,options)