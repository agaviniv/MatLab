
clear all
np=50;
roI=.5*eye(2);
Z=opmmag(1/2,'z');
ron=(Z+1/2*eye(2))/trace(Z+1/2*eye(2));
p1=ron(1,1);


Dro=zeros(2);Dro1=zeros(2); eta=zeros(1,np);
for k =1:np
    dd(k,:)=diag(roI);
    eta(k)=2*sum((dd(k,:)-.5*[1 1]).^2);
    ga=.2; %funciona como tempo!
    Dro1=atenuamp(p1,ga,'a',roI);       
    roI = Dro1;
end

subplot(2,2,1);plot(dd,'.'); axis([0 np 0 1]);
figure(2);plot(eta,'.b');


%%
clear all
np=50;
roI=.25*eye(4);
Z=opmmag(3/2,'z');
ron=(Z+3/2*eye(4))/trace(Z+3/2*eye(4));
p1=trlsq(ron);
p2=trmsq(ron);

Dro=zeros(4);Dro1=zeros(4); eta=zeros(1,np);
for k =1:np
    dd(k,:)=diag(roI);
    eta(k)=4*sum((dd(k,:)-.25*[1 1 1 1]).^2);
    ga=.2; %funciona como tempo!
    gb=.2;
    Dro1=atenuamp(p1(1,1),ga,'a',roI);
    Dro2=atenuamp(p2(1,1),gb,'b',Dro1);    
    roI = Dro2;
end

figure(1);subplot(2,2,2);plot(dd,'.');axis([0 np 0 1]);
figure(2);hold on;plot(eta,'*k');


%%
clear all
np=50;
roI=.125*eye(8);
Z=opmmag(7/2,'z');
%roI=.25*eye(8)-1e-2*Z; 
ron=(Z+3.5*eye(8))/28;
p1=trlsq(trlsq(ron));
p2=trlsq(trmsq(ron));
p3=trmsq(trmsq(ron));


Dro=zeros(8);Dro1=zeros(8); Dro2=zeros(8);
eta=zeros(1,np);
for k =1:np
    dd(k,:)=diag(roI);
    eta(k)=8*sum((dd(k,:)-.125*ones(1,8)).^2);
    ga=.2;
    gb=.2;
    gc=.2;
    Dro1=atenuamp(p1(1,1),ga,'a',roI);
    Dro2=atenuamp(p2(1,1),gb,'b',Dro1);
    Dro3=atenuamp(p3(1,1),gc,'c',Dro2);    
    roI = Dro3;
end

figure(1);subplot(2,2,3);plot(dd,'.');axis([0 np 0 1]);
figure(2);hold on; plot(eta,'+k');


%%
clear all
np=50;
roI=.0625*eye(16);
Z=opmmag(15/2,'z');
%roI=.25*eye(16)-1e-2*Z; 
ron=(Z+15/2*eye(16))/trace(Z+15/2*eye(16));
p1=trlsq(trlsq(trlsq(ron)));
p2=trlsq(trlsq(trmsq(ron))); %== trlsq(trmsq(trlsq(ron)))
p3=trmsq(trmsq(trlsq(ron))); %== trmsq(trlsq(trmsq(ron)))
p4=trmsq(trmsq(trmsq(ron)));

Dro=zeros(16);Dro1=zeros(16); Dro2=zeros(16); Dro3=zeros(16); 
eta=zeros(1,np);
for k =1:np
    dd(k,:)=diag(roI);
    eta(k)=16*sum((dd(k,:)-.0625*ones(1,16)).^2);
    ga=.2;
    gb=.2;
    gc=.2;
    gd=.2;
    Dro1=atenuamp(p1(1,1),ga,'a',roI);
    Dro2=atenuamp(p2(1,1),gb,'b',Dro1);
    Dro3=atenuamp(p3(1,1),gc,'c',Dro2);
    Dro4=atenuamp(p4(1,1),gc,'d',Dro3);
    roI = Dro3;
end

figure(1);subplot(2,2,4);plot(dd,'.');axis([0 np 0 1]);
figure(2);hold on; plot(eta,'xr'); legend('N=2','N=4','N=8','N=16')
xlabel('iteraction'); ylabel('\eta'); ylim([0 1.1])
