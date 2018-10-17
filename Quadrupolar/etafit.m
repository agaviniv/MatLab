function err = etafit(gama,t,y,p,handle)

roin = 0.25*eye(4);   
np = size(t,2);
etaf=zeros(1,np);  
for k=1:np
  if ((k ~= 2) & (k ~= 5) & (k ~= 6) & (k ~= 7)) 
    etaf(k)=trace((roin-0.25*eye(4))^2)/trace((0.25*eye(4))^2);        
  end  
    Dro = atenuamp(p(1),gama(1),'a',roin);       
    Dro2 = atenuamp(p(2),gama(2),'b',Dro);

    roin=Dro2;
end 

err = norm(etaf-y);

 set(gcf,'DoubleBuffer','on');
 set(handle,'ydata',etaf)
 drawnow
 pause(.04)
end