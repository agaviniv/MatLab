function err = disalfafit(palfa,t,y,peps,handle)

roin = 0.25*eye(4);   
np = size(t,2);
etaf=zeros(1,np);  
for k=1:np
   etaf(k)=trace((roin-0.25*eye(4))^2)/trace((0.25*eye(4))^2); 
   Dro = dissipm(peps,palfa,roin);
   roin=Dro;
end 

err = norm(etaf-y);

%  set(gcf,'DoubleBuffer','on');
%  set(handle,'ydata',wcurva)
%  drawnow
%  pause(.04)
end