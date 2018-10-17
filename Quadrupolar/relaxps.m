clear all

load DensidadePPS00W.mat;
ro00= UAB; %ro00a= UA; ro00b= UB; 
mm = max(diag(ro00{end}));
for kk=1:40; ro00{kk}=.25*eye(4)+.25*ro00{kk}/mm; pp(kk,:)=(diag(ro00{kk}));end;
figure(1);subplot(2,2,1);
plot(1e3*tempo,pp,'.'); ylim([-0.1 0.6]); %legend('\rho00','\rho01','\rho10','\rho11')

clear UA UAB UB WA WB colu fila klm rho tab1 tempo

load DensidadePPS01W.mat;
ro01= UAB; %ro01a= UA; ro01b= UB; 
mm = max(diag(ro01{end}));
for kk=1:40; ro01{kk}=.25*eye(4)+0.25*ro01{kk}/mm; pp(kk,:)=(diag(ro01{kk}));end;
subplot(2,2,2);
plot(1e3*tempo,pp,'.'); ylim([-0.1 0.6]); %legend('\rho00','\rho01','\rho10','\rho11')

clear UA UAB UB WA WB colu fila klm rho tab1 tempo

load DensidadePPS10W.mat;
ro10= UAB; %ro10a= UA; ro10b= UB; 
mm = max(diag(ro10{end}));
for kk=1:40; ro10{kk}=.25*eye(4)+0.25*ro10{kk}/mm; pp(kk,:)=(diag(ro10{kk}));end;
subplot(2,2,3);
plot(1e3*tempo,pp,'.'); ylim([-0.1 0.6]); %legend('\rho00','\rho01','\rho10','\rho11')

clear UA UAB UB WA WB colu fila klm rho tab1 tempo

load DensidadePPS11W.mat;
ro11= UAB; %ro11a= UA; ro11b= UB; 
mm = max(diag(ro11{end}));
for kk=1:40; ro11{kk}=.25*eye(4)+0.25*ro11{kk}/mm; pp(kk,:)=(diag(ro11{kk}));end;
subplot(2,2,4);
plot(1e3*tempo,pp,'.'); ylim([-0.1 0.6]); %legend('\rho00','\rho01','\rho10','\rho11')

clear UA UAB UB WA WB colu fila klm rho tab1 

%%
% ######################################################
% Calculo da contracao do espaco de fases
for kk=1:40
    roi{kk}=0.25*(ro00{kk}+ro01{kk}+ro10{kk}+ro11{kk});
    di(kk,:)=diag(roi{kk});
    eta(kk)=4*trace((roi{kk}-.25*eye(4))^2);
end
figure;
subplot(1,2,1);plot(1e3*tempo,di,'.'); ylim([-0.1 0.6]);
subplot(1,2,2);plot(1e3*tempo,eta,'.'); ylim([-0.01 0.6]);


% Fitting com canal de atenuacao de amplitude
t = 1e3*tempo'; 
figure;
plot(t,eta,'ro'); hold on; h = plot(t,eta,'b'); hold off;
title('Input data'); ylim([0 0.7])

p(1)=5/6;
p(2)=2/3;
start = [0.08;0.08]; %gama's
options = optimset('Simplex','on','TolX',0.001);
estimated_gama = fminsearch(@(x)etafit(x,t,eta,p,h),start,options)


ga = estimated_gama(1); gb = estimated_gama(2);
roi=.25*eye(4);
for k=1:40
    eta(k)=4*trace((roi - 0.25*eye(4))^2); 
    roi=atenuamp(5/6,ga,'a',roi); %
    roi=atenuamp(2/3,gb,'b',roi); %
end
hold on;
plot(tempo,eta,'r');


Slin = zeros(1,40);
for k=1:40
    Slin(k) = 1-(trace(ro00{k}^2));
end    